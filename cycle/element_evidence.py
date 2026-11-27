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
its 3' side. ISfinder entries, and the transposase ORFs of most families, run
5' to 3'; when the two disagree the ISfinder strand is used.

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


def build_isfinder_db(fasta: Path, out_dir: Path) -> Path:
    """BLAST database of the ISfinder FASTA, rebuilt when the FASTA content changes."""
    out_dir.mkdir(parents=True, exist_ok=True)
    prefix = out_dir / "isfinder"
    key_path = Path(f"{prefix}.source")
    key = hashlib.sha256(Path(fasta).read_bytes()).hexdigest()
    files = [Path(f"{prefix}{suffix}") for suffix in (".nhr", ".nin", ".nsq")]
    if all(f.exists() for f in files) and key_path.exists() and key_path.read_text() == key:
        return prefix
    key_path.unlink(missing_ok=True)
    subprocess.run(["makeblastdb", "-in", str(fasta), "-dbtype", "nucl", "-out", str(prefix)],
                   check=True, capture_output=True)
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


def annotate_sample(json_path: Path, db: Path, hmm: Path, threads: int = 1) -> Path:
    """Add ISfinder hits, transposase domains and orientation to one annotated JSON, in place."""
    records = json.loads(json_path.read_text())
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
        if by_orf and by_hit:
            status = "agree" if by_orf == by_hit else "disagree"
        else:
            status = "orf_only" if by_orf else "isfinder_only" if by_hit else "none"
        record["isfinder_hit"] = hit
        record["transposase_domains"] = sorted(names)
        record["orientation"] = by_hit or by_orf
        record["orientation_orf"] = by_orf
        record["orientation_isfinder"] = by_hit
        record["orientation_status"] = status
        if record.get("evidence") != "tail_head_reads":
            record["evidence"] = "is_reference_hit" if hit else "transposase_domain" if names else "none"

    pending = json_path.with_name(f".{json_path.name}.pending")
    pending.write_text(json.dumps(records, indent=2))
    pending.replace(json_path)
    return json_path
