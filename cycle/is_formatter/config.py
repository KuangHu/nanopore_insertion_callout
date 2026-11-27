"""Default configuration for IS element formatting / flanking extraction."""

DEFAULT_FORMATTER_OUTPUT_DIR = "data/is_formatter_output"
DEFAULT_FLANK_SIZE = 80            # bp, matches ISExtractor default
DEFAULT_MIN_READS_FOR_ASSEMBLY = 3
DEFAULT_ASSEMBLY_TIMEOUT = 120     # seconds per IS element
DEFAULT_MIN_ENTROPY = 1.7         # bits/base; filters homopolymers, tandem repeats
MAX_END_CLIP = 20                 # bases of the called consensus that may stay unaligned at each end
MAX_FALLBACK_READS = 50           # reads tried, longest first, when the local assembly gives no contig
POLISH_ROUNDS = 3                 # racon rounds on the backbone read
POLISH_PAD = 500                  # bases kept on each side of the element in the backbone
# Two records of a sample are copies of one element at this identity (by platform; the
# nanopore value allows for the errors left after polishing)...
COPY_MIN_IDENTITY = {"OXFORD_NANOPORE": 0.93, "PACBIO_SMRT": 0.95}
COPY_MIN_COVERAGE = 0.90          # ...over this fraction of both sequences
MAX_OVERHANG = 50                 # unaligned bases allowed where a polishing read leaves the backbone
