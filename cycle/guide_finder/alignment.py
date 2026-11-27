"""Short alignment finder for between-sequence comparisons.

Adapted from RNA_guide_editor_finder ShortAlignmentFinder — only the
between-sequence methods needed for guide RNA detection.
"""

from __future__ import annotations


_COMP = str.maketrans("ACGTacgt", "TGCAtgca")


def reverse_complement(seq: str) -> str:
    return seq.translate(_COMP)[::-1]


class ShortAlignmentFinder:
    """Find short DNA alignments between two sequences.

    Uses a tiered mismatch rule:
      - Exact match (0 mismatches): min_length (default 9bp)
      - With mismatches (≤ max_mismatches): min_length_for_mismatch (default 12bp)

    Parameters
    ----------
    min_length : int
        Minimum alignment length for exact (0 mismatch) hits.
    max_mismatches : int
        Maximum number of mismatches allowed for longer hits.
    min_length_for_mismatch : int
        Minimum alignment length required when mismatches > 0.
    check_forward : bool
        Find direct (forward-forward) matches.
    check_revcomp : bool
        Find reverse-complement matches.
    """

    def __init__(
        self,
        min_length: int = 9,
        max_mismatches: int = 1,
        min_length_for_mismatch: int = 12,
        check_forward: bool = True,
        check_revcomp: bool = True,
    ) -> None:
        if min_length < 1:
            raise ValueError("min_length must be at least 1")
        if max_mismatches < 0:
            raise ValueError("max_mismatches cannot be negative")
        if not check_forward and not check_revcomp:
            raise ValueError("At least one of check_forward or check_revcomp must be True")

        self.min_length = min_length
        self.max_mismatches = max_mismatches
        self.min_length_for_mismatch = min_length_for_mismatch
        self.check_forward = check_forward
        self.check_revcomp = check_revcomp

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------

    def _passes_mismatch_filter(self, length: int, mismatches: int) -> bool:
        """Check if a hit passes the tiered mismatch rule.

        - length >= min_length with 0 mismatches: always OK
        - length >= min_length_for_mismatch with <= max_mismatches: OK
        - otherwise: rejected
        """
        if mismatches == 0 and length >= self.min_length:
            return True
        if mismatches <= self.max_mismatches and length >= self.min_length_for_mismatch:
            return True
        return False

    def find_alignments_between(
        self, sequence1: str, sequence2: str
    ) -> list[dict]:
        """Find all short alignments between two different sequences.

        Parameters
        ----------
        sequence1 : str
            First DNA sequence (e.g., flanking region).
        sequence2 : str
            Second DNA sequence (e.g., noncoding region).

        Returns
        -------
        list[dict]
            Each dict contains: pos1, pos2, length, seq1, seq2,
            mismatches, mismatch_positions, orientation.
        """
        sequence1 = sequence1.upper()
        sequence2 = sequence2.upper()

        found = []
        if self.check_forward:
            found += self._scan_diagonals(sequence1, sequence2, "forward")
        if self.check_revcomp:
            # Compare with the reverse complement, then give positions on sequence2 itself
            for hit in self._scan_diagonals(sequence1, reverse_complement(sequence2), "reverse_complement"):
                hit["pos2"] = len(sequence2) - hit["pos2"] - hit["length"]
                hit["seq2"] = sequence2[hit["pos2"]:hit["pos2"] + hit["length"]]
                found.append(hit)

        consolidated = self._consolidate_matches_between(found)
        consolidated.sort(key=lambda x: (x["pos1"], x["pos2"]))
        return consolidated

    # ------------------------------------------------------------------
    # Diagonal scan
    # ------------------------------------------------------------------

    def _scan_diagonals(self, sequence1: str, sequence2: str, orientation: str) -> list[dict]:
        """All ungapped alignments of the two sequences that pass the tiered rule.

        Every offset between the sequences is one diagonal. On a diagonal the
        mismatch positions cut it into exact runs; a window may span up to
        max_mismatches of them. Each longest window that starts and ends on a
        match is reported if it passes the rule, and so is each exact run that
        is long enough by itself, because a run can pass where the longer
        window around it, with its mismatch, is too short for the mismatch tier.
        """
        hits = []
        len1, len2 = len(sequence1), len(sequence2)
        for offset in range(-(len1 - 1), len2):
            start1 = max(0, -offset)
            start2 = start1 + offset
            size = min(len1 - start1, len2 - start2)
            if size < self.min_length:
                continue
            mismatches = [i for i in range(size) if sequence1[start1 + i] != sequence2[start2 + i]]
            bounds = [-1] + mismatches + [size]
            windows = set()
            for i in range(len(bounds) - 1):
                windows.add((bounds[i] + 1, bounds[i + 1]))  # exact run
                windows.add((bounds[i] + 1, bounds[min(i + 1 + self.max_mismatches, len(bounds) - 1)]))
            for begin, end in windows:
                # Trim mismatches from the ends (two neighbouring mismatches leave one there)
                while begin < end and sequence1[start1 + begin] != sequence2[start2 + begin]:
                    begin += 1
                while end > begin and sequence1[start1 + end - 1] != sequence2[start2 + end - 1]:
                    end -= 1
                inside = [m - begin for m in mismatches if begin <= m < end]
                if not self._passes_mismatch_filter(end - begin, len(inside)):
                    continue
                hits.append({
                    "pos1": start1 + begin,
                    "pos2": start2 + begin,
                    "length": end - begin,
                    "seq1": sequence1[start1 + begin:start1 + end],
                    "seq2": sequence2[start2 + begin:start2 + end],
                    "mismatches": len(inside),
                    "mismatch_positions": inside,
                    "orientation": orientation,
                })
        return hits

    # ------------------------------------------------------------------
    # Consolidation
    # ------------------------------------------------------------------

    def _consolidate_matches_between(self, matches: list[dict]) -> list[dict]:
        if not matches:
            return []

        matches_sorted = sorted(
            matches, key=lambda x: (-x["length"], x["mismatches"])
        )
        consolidated: list[dict] = []

        for match in matches_sorted:
            if not any(
                self._matches_overlap_between(match, ex) for ex in consolidated
            ):
                consolidated.append(match)

        return consolidated

    def _matches_overlap_between(self, m1: dict, m2: dict) -> bool:
        if m1["orientation"] != m2["orientation"]:
            return False
        s1_ovl = self._regions_overlap(
            m1["pos1"], m1["pos1"] + m1["length"],
            m2["pos1"], m2["pos1"] + m2["length"],
        )
        s2_ovl = self._regions_overlap(
            m1["pos2"], m1["pos2"] + m1["length"],
            m2["pos2"], m2["pos2"] + m2["length"],
        )
        return s1_ovl and s2_ovl

    # ------------------------------------------------------------------
    # Utilities
    # ------------------------------------------------------------------

    @staticmethod
    def _regions_overlap(s1: int, e1: int, s2: int, e2: int) -> bool:
        return not (e1 <= s2 or e2 <= s1)
