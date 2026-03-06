"""Taxonomic profiling with MetaPhlAn4.

Generates species-level relative abundance profiles from metagenomic reads
using clade-specific marker genes.

TODO:
    - [ ] Implement MetaPhlAn4 wrapper
    - [ ] Parse output into standardized abundance matrix
    - [ ] Add merged abundance table across samples
"""

import logging

logger = logging.getLogger(__name__)


def run_metaphlan(
    fastq: str,
    output_dir: str,
    database: str = "databases/metaphlan4",
    threads: int = 8,
) -> dict:
    """Run MetaPhlAn4 taxonomic profiling on a sample.

    Args:
        fastq: Path to input FASTQ file.
        output_dir: Directory for abundance profiles.
        database: Path to MetaPhlAn4 Bowtie2 database.
        threads: Number of threads.

    Returns:
        Dictionary with abundance profile path and species count.
    """
    raise NotImplementedError(
        "MetaPhlAn4 integration is planned for Phase 3. "
        "See README.md roadmap for current progress."
    )
