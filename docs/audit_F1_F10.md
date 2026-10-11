# Audit findings F1–F10 — insertion / flanking correctness

**Baseline commit audited:** `90719c6ace42bf7f342ae3ab672247bc8a14bc7c`
**Audit scope:** Sniffles VCF → insertion table → circle bait / reads → formatter read-recruitment, assembly, locate, fallback paths, JSON/TSV outputs; CLI defaults, batch success semantics, downstream coordinate consumers.
**What this doc is for:** a durable record of the correctness gaps the external review identified, their current remediation status, and the planned fix order. The audit did not modify business code; neither does this doc. Any analysis using current outputs must consult the "Impact on existing deliverables" table before being treated as validated.

---

## Severity + status overview

| ID  | Severity | Area                     | Status                          |
|-----|----------|--------------------------|---------------------------------|
| F1  | P1       | sniffles_runner `<INS>`  | open                            |
| F2  | P1       | read recruitment anchor  | open                            |
| F3  | P1       | formatter locate quality | partially addressed (v3 refiner only; formatter itself unchanged) |
| F4  | P1       | TH-junction counting     | open                            |
| F5  | P2       | coordinate semantics     | partially addressed (one downstream consumer patched) |
| F6  | P2       | strand propagation       | open                            |
| F7  | P2       | success/failure labeling | partially addressed (v3 refiner records status; formatter unchanged) |
| F8  | P2       | rerun cache consistency  | open                            |
| F9  | P2       | candidate fallback       | open                            |
| F10 | P2       | assembly quality label   | open                            |

Status key:
* **open** — not touched in the current repo; outputs cannot be trusted for this concern.
* **partially addressed** — a specific downstream consumer or newer catalog branch mitigates the symptom, but the upstream source of the error is unchanged. Not safe to extrapolate the mitigation to callers that still read raw formatter output.

---

## Findings

### F1 [P1] Symbolic `<INS>` emits entire assembly contig as insertion
**Location:** `cycle/preprocess/sniffles_runner.py:142-154`, `270-355`.
**Behaviour:** For `<INS>` SV records (no explicit ALT, SVLEN only), `_assemble_insertion` extracts supporting reads, assembles, returns the longest unitig unchanged. `run_sample` writes that unitig into the record's `Sequence` field. No re-anchoring to flanking reference, no clipping of the inserted segment.
**Consequence:** Insertion consensus can contain substantial host DNA on either side of the real inserted segment. All downstream HMM / 6-frame / IR / guide analyses on those records are polluted. `[IS][IS]` bait built from this consensus is also polluted.
**Reproduction (from audit):** With declared SVLEN=1000 and a simulated local assembly returning 10 kb, the final `Sequence` is written as the full 10 kb.
**Planned fix:** Use reference-anchor + supporting-read insertion CIGAR to localize the inserted interval, clip to it, validate. Records whose inserted segment cannot be isolated must be kept as `unresolved`, not labelled as insertion.
**Status:** open.

### F2 [P1] Read recruitment lacks site identity — flanks can come from a sister copy or circle
**Location:** `cycle/circle_detect/circle_finder.py:135-172`, `206-239`; `cycle/is_formatter/formatter.py:149-204`, `585-595`.
**Behaviour:** Bait is just `insertion_seq ++ insertion_seq` with no site-specific genomic anchor. Formatter groups all primary reads mapping to the bait by UUID; no MAPQ threshold, no original-genome-position check, no requirement that reads be the Sniffles-named supporting reads.
**Consequence:** When the same IS element exists at multiple genomic loci (routine for IS110), reads from locus A can be recruited into the record for locus B. Flanks written for locus B are then from locus A — site identity is not established.
**Reproduction (from audit):** MAPQ=0 bait reads enter assembly input; a candidate contig of form `IS+IS+IS` can have its middle copy located and the outer IS pieces returned as "flanks".
**Planned fix:** Build site identity from the original-genome BAM (require ≥1 uniquely-anchored supporting read per locus). Separate circle (episomal) molecules from chromosomal integrants. Mark homologous-copy ambiguity as `ambiguous`, not as a validated locus flank.
**Status:** open.

### F3 [P1] Local or supplementary alignments determine full insertion boundaries
**Location:** `cycle/is_formatter/formatter.py:365-396`.
**Behaviour:** The locate step filters only `flag & 4`; picks the single highest-MAPQ record. No query coverage check, no end-clipping check, no identity check, no minimum MAPQ, no verification that split alignments are complete. Supplementary records compete on equal footing with primary.
**Consequence:** A 1 kb query with CIGAR `400S200M400S` is called a 200 bp insertion with both ends truncated. MAPQ=0 accepted. A high-MAPQ `400H200M400H` supplementary can override a complete primary and define a wrong boundary. ORF/guide/junction-motif downstream all inherit the error.
**Mitigation present:** The v3 junction refiner (`junction_refiner_v3.py`, `assembly_refiner_v3.py`) does enforce coverage/identity/clipping via Type A (single-alignment inline I ≥50 bp) + Type B (primary + SA split, both from the same read) + MAPQ ≥20 filter, and emits `junction_precision` + `junction_architecture` labels. This mitigates downstream consumers that use v3 records. The formatter's own `locate()` function is unchanged.
**Planned fix:** Apply the same thresholds inside `formatter.py:locate()` so the primary output is correct, not only the v3 catalog.
**Status:** partially addressed.

### F4 [P1] Deletions crossing the junction counted as tail-head evidence
**Location:** `cycle/circle_detect/circle_finder.py:343-346`.
**Behaviour:** TH call fires whenever `reference_start < N < reference_end`. Reference span through a D/N CIGAR op does not mean the read's query bases actually bridge the junction.
**Consequence:** The TH count that drives `require_th_reads` is false-positive-prone. The default `require_th_reads=True` activity filter accepts candidates that have no real junction read support.
**Reproduction (from audit):** With N=1000, `reference_start=800`, CIGAR=`100M200D100M`, `reference_end=1200`, this is counted as 1 TH read even though the junction sits inside the deletion.
**Impact on circle_evidence counts:** likely over-estimated — magnitude unknown.
**Planned fix:** Walk CIGAR near the junction, require real match anchors on both sides of N, reject large indels straddling the junction, require query-base coverage on both sides. Add this specific negative case to regression tests.
**Status:** open.

### F5 [P2] Insertion length written as a reference interval
**Location:** `cycle/preprocess/sniffles_runner.py:253-256`; `cycle/is_formatter/formatter.py:462-478`.
**Behaviour:** Sniffles table sets `End = POS + SVLEN`; this is copied verbatim into `is_element.end` and `insertion_site.pos_3p`. For pure INS, SVLEN is the length of the newly inserted sequence — it is not a reference coordinate. Adding it to POS produces a non-existent right breakpoint on the reference.
**Reproduction:** POS=1000, SVLEN=1000 INS → output `chr1:1000-2000`, `pos_3p=2000`. That metadata does not represent any real reference breakpoint.
**Mitigation present:** The flanking extractor in `scripts/build_all_insertions_with_flanking.py` (lives in the sibling `cycle_rna_analysis` tree, not in this repo) was patched to use `down_s = ins_start` instead of `down_s = ins_end`, with the fix locked via a `flank_extractor_version = v2_ins_point_coordinate` tag. That single consumer now reads downstream flanks correctly. `is_element.end` and `insertion_site.pos_3p` elsewhere in the pipeline are still the wrong coordinate.
**Planned fix:** Separate reference breakpoint (coordinate system + uncertainty) from assembly-internal interval from insertion length in the schema. Model DUP/BND distinctly from INS.
**Status:** partially addressed.

### F6 [P2] Relative strand dropped; output fixed to `+`
**Location:** `cycle/is_formatter/formatter.py:437-442`, `462-465`, `590-595`.
**Behaviour:** `locate` returns strand, but neither the assembly path nor the longest-read path persists it. Final record always has `strand = +`.
**Clarification (from audit):** This does not necessarily break `upstream + IS + downstream` concatenation and does not necessarily invalidate ORFs re-predicted on the output sequence. The problem is: downstream code cannot assume the sequence orientation matches the original consensus or the reference.
**Planned fix:** Define an explicit output orientation convention. If normalizing to a consensus direction, reverse-complement the whole record and swap flanks when needed; also store orientation relative to the reference. Do not mix the two.
**Status:** open.

### F7 [P2] Failures / partial results indistinguishable from validated results in output
**Location:** `cycle/is_formatter/formatter.py:398-412`, `590-622`, `689-695`.
**Behaviour:** Both extraction paths can fail yet still `append` a record and return the directory. Batch-success logic treats "non-empty directory" as success. Flank slicing allows zero length. Assembly path tags `assembly_method=assembly` whenever an alignment exists, regardless of whether both ends resolved.
**Reproduction:** With no reads: empty insertion, `assembly_method=none`. With a contig containing only IS: both flanks at 0 bp but `assembly_method=assembly`.
**Current partial help:** `assembly_method` and the two flank-length fields let a careful human spot some failures. But sample-level success status and the normal record stream are not separated — downstream code cannot filter on validity without parsing field heuristics.
**Mitigation present:** The v3 catalogs (`sniffles_catalog_v3`, `assembly_catalog_v3`) do emit explicit `status` values (`refined`, `insufficient_reads`, `ref_missing`, `fetch_fail`). Formatter's primary output is unchanged.
**Planned fix:** Per-record `status ∈ {failed, partial, validated, ambiguous}`; keep failure diagnostics but segregate them from trusted records; emit per-side flank length and support count; batch report with success/incomplete/failed tallies.
**Status:** partially addressed.

### F8 [P2] Rerun inconsistency — stale bait/BAM mixed with new candidates
**Location:** `cycle/circle_detect/circle_finder.py:151-153`, `217-219`; `cycle/preprocess/sniffles_runner.py:253-254`.
**Behaviour:** Bait and BAM are reused based solely on file existence. Re-parsing the VCF generates fresh random UUIDs. If someone rebuilds the insertion table but keeps circle outputs, old bait UUIDs no longer match new entries — junction statistics will silently skip the old IDs and can go to zero.
**Reproduction:** An existing `old__th__j1000` bait keeps returning the `old` mapping even after being passed a new ID.
**Planned fix:** Cache keys must bind to input content + candidate ID + parameter fingerprint. Rebuild dependencies consistently on upstream change. Stable IDs alone do not solve consensus-content changes.
**Status:** open.

### F9 [P2] Only tries longest read / contig; gives up when it fails to locate
**Location:** `cycle/is_formatter/formatter.py:306-327`, `429-435`.
**Behaviour:** Candidates ranked by total length, not by target coverage or dual-flank support. If the longest read fails to locate IS, no retry on the next read. If the longest contig does not contain the full IS, no search of other unitigs. Partial-locate success prevents fallback.
**Reproduction:** Longest read locate fails, second-longest read would locate cleanly, pipeline still returns `None` with only one locate attempt.
**Planned fix:** Rank candidates by full-target coverage + dual anchors + locus uniqueness before picking. Record rejection reasons for the losing candidates.
**Status:** open.

### F10 [P2] Polished vs raw assembly not distinguished in quality label
**Location:** `cycle/is_formatter/formatter.py:278-304`, `595`; `cycle/preprocess/sniffles_runner.py:327-355`.
**Behaviour:** When minipolish is missing or fails, the raw miniasm unitig is returned and the formatter still tags the record `assembly_method=assembly`. Symbolic-INS upstream assembly path also doesn't polish.
**Why this matters:** miniasm does not produce consensus — unitig base error rate is close to the input reads' error rate (~5-10% for nanopore). An `assembly` tag therefore does not imply consensus-grade sequence. ORF/short-motif analyses that assume consensus accuracy are invalid on raw unitigs.
**Planned fix:** Distinguish `single_read` / `raw_assembly` / `polished_assembly`; record polishing status, supporting-read count, and verification metrics. Short-motif / ORF analyses restricted to records meeting a quality bar.
**Status:** open.

---

## Confirmed correct (not findings; do not re-flag)

1. SAM `POS` → 0-based conversion by subtracting 1; reference-end computed from ref-consuming CIGAR operations.
2. Flank slicing is off-by-one-clean given correct boundaries and strand.
3. Formatter deduplicates by read name within the same UUID and excludes secondary/supplementary records; helps avoid double-counting but does not establish site identity.
4. BAM `query_sequence` is stored in alignment orientation; absence of `get_forward_sequence` is not a universal reverse-complement bug here.
5. `require_th_reads=True` is an explicit task filter, not a general insertion-discovery completeness policy. Changing it only changes candidate scope; it does not fix any of the quality issues above.
6. Prior fix for `orf_annotation=None` (commit `90719c6`) is valid and independent of the audit concerns.

---

## Impact on existing deliverables

| Deliverable | Findings touched | Current trust level |
|---|---|---|
| `active_is110_000_094` (2,475 IS110 records) | F1, F2, F4, F6, F10 | Candidate list — do not treat flanks or coordinates as validated |
| `active_is110_006_078` (1,926 IS110, prior snapshot) | same | same |
| `is110_consensus_000_084_6frame.fna` (18,590 HMM-matched IS110) | F1 (if any consensus is contaminated), F2 | likely sufficient for HMM classification; not safe for short-motif / junction biology |
| `/groups/.../all_insertions_fixed_000_094/` (Sniffles v2 shards, `v2_ins_point_coordinate`) | F5 fix applied for downstream flank only | downstream flank correct for extraction; upstream pipeline still has F1, F2, F3, F6, F10 |
| `all_insertions_refined_v3/` (currently running) | F3 mitigated per record via Type A+B refinement; F1/F2 still inherited | the exact-anchor + sequence-identity ≥0.9 subset (~15K) is the safest pool |
| `assembly_catalog_v3/` (currently queued) | same mitigation as above | inherits upstream F1/F2 |
| LTG/RTG, cargo, endpoint-sharing analyses | all upstream findings propagate | do not re-run against current outputs until F1/F2/F3 patched |

---

## Fix order and acceptance

**Phase 1 — stop trusting flanks until these pass**: F1, F2, F3, F4.
**Phase 2 — status + semantics hygiene**: F5, F7, F10.
**Phase 3 — robustness**: F6, F8, F9.

**Acceptance must cover at minimum:**
- Known insertions in both orientations
- Dual 80 bp flanks and single-sided misses
- Same IS at multiple loci with distinct flanks
- Tandem / circular reads
- Partial matches + supplementary alignments
- Large deletion at junction
- `<INS>` reconstruction
- Assembly / polishing failures
- Rerun after candidate changes

**Metrics:** on synthetic / known-truth data report insertion-boundary offset, flank-locus accuracy, base identity, complete-output fraction, and `ambiguous`/`failed` fractions. The audit's 15 diagnostic checks should be converted into expected-behavior assertions after fixes land; current all-OK is not acceptance.

---

## References

- SAM spec: https://github.com/samtools/hts-specs/blob/master/SAMv1.tex
- VCF spec: https://github.com/samtools/hts-specs/blob/master/VCFv4.3.tex
- miniasm: https://github.com/lh3/miniasm
- Audit reproduction script lives at `review_checks/audit_flanks.py` in the review repo; standard library only, does not touch business source.
