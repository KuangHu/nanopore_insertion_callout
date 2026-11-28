#!/usr/bin/env python3
"""Move IS110 element intervals between the 3' and 5' placements of their core.

An IS110 element shares a short core (usually 2 bp) with its target and can be
placed k bases to either side with the same sequence. Cycle records follow
ISfinder and put the shared bases at the element's 5' end; Cross_reference_IS
reports the rightmost placement in the element's direction, which puts them at
the 3' end.

Usage:
    python scripts/convert_core_convention.py --genome filled.fna --to 5prime in.tsv > out.tsv

in.tsv needs the columns contig, start, end (0-based, half-open interval of the
element in --genome, the genome that holds it) and strand (+ or -, the
element's direction). Every column is kept; start_<to>, end_<to> and core_shift
(reference bases, negative towards the contig start) are added. Runs of more
than MAX_CORE shared bases are left where they are (shift 0).
"""

from __future__ import annotations

import argparse
import csv
import sys

MAX_CORE = 10  # as element_evidence.MAX_END_REPEAT


def read_fasta(path: str) -> dict[str, str]:
    sequences, name, parts = {}, None, []
    with open(path) as handle:
        for line in handle:
            if line.startswith(">"):
                if name is not None:
                    sequences[name] = "".join(parts).upper()
                name, parts = line[1:].split()[0], []
            else:
                parts.append(line.strip())
    if name is not None:
        sequences[name] = "".join(parts).upper()
    return sequences


def core_shift(contig: str, start: int, end: int, strand: str, to: str) -> int:
    """Reference bases to move the interval so the shared bases end up at the element's ``to`` end.

    Bases that the element's reference-right end shares with the bases before it move the
    element left; bases its reference-left end shares with the bases after it move it right.
    In the element's direction the reference-right end is its 3' end on the plus strand.
    """
    move_left = (strand == "+") == (to == "5prime")
    k = 0
    if move_left:
        while k < end - start and start - 1 - k >= 0 and contig[start - 1 - k] == contig[end - 1 - k]:
            k += 1
    else:
        while k < end - start and end + k < len(contig) and contig[end + k] == contig[start + k]:
            k += 1
    if k > MAX_CORE:
        return 0
    return -k if move_left else k


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("table")
    parser.add_argument("--genome", required=True)
    parser.add_argument("--to", choices=["5prime", "3prime"], default="5prime")
    args = parser.parse_args()

    genome = read_fasta(args.genome)
    with open(args.table, newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields = reader.fieldnames + [f"start_{args.to}", f"end_{args.to}", "core_shift"]
        writer = csv.DictWriter(sys.stdout, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in reader:
            start, end = int(row["start"]), int(row["end"])
            shift = core_shift(genome[row["contig"]], start, end, row["strand"], args.to)
            row.update({f"start_{args.to}": start + shift, f"end_{args.to}": end + shift, "core_shift": shift})
            writer.writerow(row)


if __name__ == "__main__":
    main()
