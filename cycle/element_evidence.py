"""Reference and protein evidence for each element, and its orientation.

Runs on ``*_is_records_annotated.json`` (after the ORF annotator) and adds to
every record:

  isfinder_hit          best blastn hit against ISfinder (task blastn, not
                        megablast), or None
  transposase_domains   Pfam transposase domains found on the record's ORFs
  orientation           strand of the element's 5' to 3' direction on the
                        reference: "+", "-" or "" when unknown
  orientation_orf       strand of the ORF with the best transposase domain
  orientation_isfinder  strand of the ISfinder hit
  orientation_status    agree, disagree, orf_only, isfinder_only or none

Record sequences are on the reference strand, so a minus orientation means
the element's 5' end is at the record's end and ``flanking_upstream`` lies on
its 3' side. Non-coding regions get ``element_end`` in the element's own
direction (the ORF annotator's ``type`` follows the record). Orientation
follows the ISfinder entry, which is stored from the element's left end.

Some entries carry their transposase on the minus strand (many IS200/IS605,
some IS3, IS256, IS110), so the ORF strand is compared with the strand the
entry's own transposase ORF predicts: found once per database with prodigal
and the same Pfam models. ``orientation_status`` is "two_elements" when the
best transposase ORF lies outside the ISfinder hit: the record then holds
another element beside the matched one and has no single orientation (the
ISfinder strand is kept).

An IS110 element (ISfinder family IS110, or an IS110 transposase domain
without a hit) inserts without a target duplication but shares a 2 bp core
with the target, so it has several placements with the same sequence. Such
records are moved to the placement that puts the shared bases at the
element's 5' end, as ISfinder entries store them. Other families are left
alone: there the shared bases are the target duplication;
``end_repeat_shift`` gives the move in reference bases (negative: towards the
reference start), and ``boundary_convention`` names the placement used:
"is110_core_5prime" for IS110 records, "as_called" for the rest. Cross_reference_IS
puts the core at the 3' end instead; scripts/convert_core_convention.py converts
coordinates between the two. Coordinates inside the record (ORFs, non-coding regions,
the ISfinder hit) move with it.

``evidence`` stays "tail_head_reads" when the record has them; otherwise it
becomes "is_reference_hit" for an ISfinder hit, "transposase_domain" for a
domain only, or stays "none".
"""

from __future__ import annotations

import hashlib
import json
import logging
import subprocess
import tempfile
from pathlib import Path

logger = logging.getLogger(__name__)

MIN_HIT_COVERAGE = 0.5   # fraction of the record covered by the best ISfinder hit
MAX_HIT_EVALUE = 1e-10
MAX_END_REPEAT = 10      # longer runs shared with the flank stay where they were called
IS110_DOMAINS = {"DEDD_Tnp_IS110", "Transposase_20"}


def build_isfinder_db(fasta: Path, out_dir: Path, hmm: Path, threads: int = 1) -> Path:
    """BLAST database of the ISfinder FASTA and each entry's transposase strand.

    Rebuilt when the FASTA or the HMM file changes. The strands go to
    ``<prefix>.tnp_strand.tsv``: entry name, strand of its best-scoring
    transposase ORF.
    """
    out_dir.mkdir(parents=True, exist_ok=True)
    prefix = out_dir / "isfinder"
    key_path = Path(f"{prefix}.source")
    key = hashlib.sha256(Path(fasta).read_bytes() + Path(hmm).read_bytes()).hexdigest()
    files = [Path(f"{prefix}{suffix}") for suffix in (".nhr", ".nin", ".nsq", ".tnp_strand.tsv")]
    if all(f.exists() for f in files) and key_path.exists() and key_path.read_text() == key:
        return prefix
    key_path.unlink(missing_ok=True)
    subprocess.run(["makeblastdb", "-in", str(fasta), "-dbtype", "nucl", "-out", str(prefix)],
                   check=True, capture_output=True)
    with tempfile.TemporaryDirectory() as tmp:
        work = Path(tmp)
        subprocess.run(["prodigal", "-p", "meta", "-q", "-i", str(fasta), "-a", str(work / "orfs.faa"),
                        "-o", "/dev/null"], check=True, capture_output=True)
        proteins, strands, name = {}, {}, None
        for line in (work / "orfs.faa").read_text().splitlines():
            if line.startswith(">"):
                fields = line[1:].split(" # ")
                name = fields[0].split()[0]
                strands[name] = "+" if fields[3] == "1" else "-"
                proteins[name] = ""
            else:
                proteins[name] += line
        best: dict[str, tuple[float, str]] = {}
        for orf, found in _transposase_domains(proteins, hmm, work, threads).items():
            entry = orf.rsplit("_", 1)[0]
            score = max(s for _, s in found)
            if entry not in best or score > best[entry][0]:
                best[entry] = (score, strands[orf])
    Path(f"{prefix}.tnp_strand.tsv").write_text("".join(f"{e}\t{strand}\n" for e, (_, strand) in sorted(best.items())))
    key_path.write_text(key)
    return prefix


def _isfinder_hits(sequences: dict[str, str], db: Path, work: Path, threads: int) -> dict[str, dict]:
    """Best hit per record by bit score, kept when it covers MIN_HIT_COVERAGE of the record."""
    query = work / "records.fna"
    query.write_text("".join(f">{name}\n{seq}\n" for name, seq in sequences.items()))
    result = subprocess.run(
        ["blastn", "-task", "blastn", "-query", str(query), "-db", str(db),
         "-outfmt", "6 qseqid sseqid pident length qstart qend sstart send evalue bitscore slen",
         "-evalue", str(MAX_HIT_EVALUE), "-max_target_seqs", "5", "-num_threads", str(threads)],
        check=True, capture_output=True, text=True,
    )
    best: dict[str, dict] = {}
    for line in result.stdout.splitlines():
        q, s, pident, length, qs, qe, ss, se, evalue, bits, slen = line.split("\t")
        coverage = (int(qe) - int(qs) + 1) / len(sequences[q])
        if coverage < MIN_HIT_COVERAGE:
            continue
        if q in best and float(bits) <= best[q]["bitscore"]:
            continue
        parts = s.split("_")
        best[q] = {
            "name": s,
            "query_start": int(qs),
            "query_end": int(qe),
            "family": parts[-2] if len(parts) >= 3 else "",
            "group": parts[-1] if len(parts) >= 3 else "",
            "pident": float(pident),
            "query_coverage": round(coverage, 4),
            "subject_coverage": round((abs(int(se) - int(ss)) + 1) / int(slen), 4),
            "strand": "+" if int(ss) < int(se) else "-",
            "evalue": float(evalue),
            "bitscore": float(bits),
        }
    return best


def _transposase_domains(proteins: dict[str, str], hmm: Path, work: Path, threads: int) -> dict[str, list]:
    """Pfam transposase domains per protein, with gathering thresholds: name -> [(domain, score)]."""
    fasta = work / "proteins.faa"
    fasta.write_text("".join(f">{name}\n{seq.rstrip('*')}\n" for name, seq in proteins.items()))
    table = work / "domains.tbl"
    subprocess.run(["hmmsearch", "--cut_ga", "--cpu", str(threads), "--tblout", str(table),
                    "-o", "/dev/null", str(hmm), str(fasta)], check=True, capture_output=True)
    found: dict[str, list] = {}
    for line in table.read_text().splitlines():
        if line.startswith("#"):
            continue
        fields = line.split()
        found.setdefault(fields[0], []).append((fields[2], float(fields[5])))
    return found


def _place_end_repeat(record: dict) -> int:
    """Move the element so bases it shares with a flank sit at its 5' end; returns the shift.

    A plus element whose last k bases equal the last k bases of the upstream flank can start k
    bases earlier with the same sequence; a minus element whose first k bases equal the first k of
    the downstream flank can start k bases later. Both put the shared bases at the 5' end.
    """
    element = record["is_element"]["sequence"]
    up = (record.get("flanking_upstream") or {}).get("sequence", "")
    down = (record.get("flanking_downstream") or {}).get("sequence", "")
    k = 0
    if record["orientation"] == "+":
        while k < min(len(up), len(element)) and up[-1 - k] == element[-1 - k]:
            k += 1
    elif record["orientation"] == "-":
        while k < min(len(down), len(element)) and down[k] == element[k]:
            k += 1
    if not 0 < k <= MAX_END_REPEAT:
        return 0
    if record["orientation"] == "+":
        shift, element, up, down = -k, up[-k:] + element[:-k], up[:-k], element[-k:] + down
    else:
        shift, element, up, down = k, element[k:] + down[:k], up + element[:k], down[k:]
    record["is_element"]["sequence"] = element
    for key in ("start", "end"):
        record["is_element"][key] += shift
    site = record.get("insertion_site") or {}
    for key in ("pos_5p", "pos_3p"):
        if key in site:
            site[key] += shift
    record["flanking_upstream"].update(sequence=up, length=len(up))
    record["flanking_downstream"].update(sequence=down, length=len(down))

    # 1-based positions in the record move the other way; spans reaching an end stay at that end
    n = len(element)
    def move(start, end):
        return (1 if start == 1 else min(max(start - shift, 1), n), n if end == n else min(max(end - shift, 1), n))
    annotation = record.get("orf_annotation") or {}
    for span in annotation.get("orfs", []) + (annotation.get("noncoding_regions") or []):
        span["start"], span["end"] = move(span["start"], span["end"])
        if "length" in span:
            span["length"] = span["end"] - span["start"] + 1
    if record.get("isfinder_hit"):
        record["isfinder_hit"]["query_start"], record["isfinder_hit"]["query_end"] = move(
            record["isfinder_hit"]["query_start"], record["isfinder_hit"]["query_end"])
    return shift


def annotate_sample(json_path: Path, db: Path, hmm: Path, threads: int = 1) -> Path:
    """Add ISfinder hits, transposase domains and orientation to one annotated JSON, in place."""
    records = json.loads(json_path.read_text())
    entry_strand = {}
    strand_table = Path(f"{db}.tnp_strand.tsv")
    if strand_table.exists():
        entry_strand = dict(line.split("\t") for line in strand_table.read_text().splitlines())
    sequences = {str(i): r["is_element"]["sequence"] for i, r in enumerate(records) if r["is_element"]["sequence"]}
    proteins = {}
    for i, record in enumerate(records):
        for j, orf in enumerate((record.get("orf_annotation") or {}).get("orfs", [])):
            if orf.get("protein_sequence"):
                proteins[f"{i}_{j}"] = orf["protein_sequence"]

    with tempfile.TemporaryDirectory() as tmp:
        work = Path(tmp)
        hits = _isfinder_hits(sequences, db, work, threads) if sequences else {}
        domains = _transposase_domains(proteins, hmm, work, threads) if proteins else {}

    for i, record in enumerate(records):
        hit = hits.get(str(i))
        best_orf, best_score, names = None, 0.0, set()
        for j, orf in enumerate((record.get("orf_annotation") or {}).get("orfs", [])):
            found = domains.get(f"{i}_{j}", [])
            orf["transposase_domains"] = sorted({name for name, _ in found})
            names.update(orf["transposase_domains"])
            for _, score in found:
                if score > best_score:
                    best_orf, best_score = orf, score
        by_orf = best_orf["strand"] if best_orf else ""
        by_hit = hit["strand"] if hit else ""
        if by_orf and by_hit and (best_orf["end"] < hit["query_start"] or best_orf["start"] > hit["query_end"]):
            status = "two_elements"
        elif by_orf and by_hit:
            # The ORF strand the hit predicts: the hit strand, turned over when the entry's own
            # transposase is on its minus strand.
            expected = by_hit if entry_strand.get(hit["name"], "+") == "+" else {"+": "-", "-": "+"}[by_hit]
            status = "agree" if by_orf == expected else "disagree"
        else:
            status = "orf_only" if by_orf else "isfinder_only" if by_hit else "none"
        record["isfinder_hit"] = hit
        record["transposase_domains"] = sorted(names)
        record["orientation"] = by_hit or by_orf
        record["orientation_orf"] = by_orf
        record["orientation_isfinder"] = by_hit
        record["orientation_status"] = status
        is110 = (hit or {}).get("family") == "IS110" or (not hit and names & IS110_DOMAINS)
        record["end_repeat_shift"] = _place_end_repeat(record) if is110 and record["is_element"]["sequence"] else 0
        record["boundary_convention"] = "is110_core_5prime" if is110 else "as_called"
        # The ORF annotator names non-coding regions by the record's own direction; give the
        # element's ends as well. A minus element starts at the record's end.
        swap = {"5_prime": "3_prime", "3_prime": "5_prime"}
        for region in (record.get("orf_annotation") or {}).get("noncoding_regions") or []:
            kind = region.get("type", "")
            region["element_end"] = "" if not record["orientation"] else swap.get(kind, kind) if record["orientation"] == "-" else kind
        if record.get("evidence") != "tail_head_reads":
            record["evidence"] = "is_reference_hit" if hit else "transposase_domain" if names else "none"

    pending = json_path.with_name(f".{json_path.name}.pending")
    pending.write_text(json.dumps(records, indent=2))
    pending.replace(json_path)
    return json_path
