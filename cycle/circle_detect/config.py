"""Default configuration for circular intermediate detection."""

DEFAULT_CIRCLE_OUTPUT_DIR = "data/circle_output"
DEFAULT_MIN_JUNCTION_OVERLAP = 100  # bp on each side of junction
DEFAULT_MIN_CONSENSUS_LENGTH = 200  # skip IS elements shorter than this
DEFAULT_BOUNDARY_TOLERANCE = 50  # bp tolerance for IS copy boundary detection
DEFAULT_MIN_CONSENSUS_ENTROPY = 1.7  # bits/base; filters homopolymers, tandem repeats, all-Ns
DEFAULT_MAX_TARGET_DUPLICATION = 20  # longest junction deletion accepted as a target site duplication
DEFAULT_MIN_DUPLICATION_READS = 3  # tail-head reads that must agree before a target site duplication is recorded
DEFAULT_JUNCTION_WINDOW = 5  # deletions this close to the bait junction count towards the duplication
