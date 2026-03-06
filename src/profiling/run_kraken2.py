"""Taxonomic profiling with Kraken2 and Bracken abundance estimation.

Classifies metagenomic reads against a reference database to determine
which microbial species are present and estimate their relative abundance.

Example:
    >>> from src.profiling.run_kraken2 import classify_reads, estimate_abundance
    >>> report = classify_reads("sample.fastq.gz", "databases/kraken2_standard/")
    >>> abundances = estimate_abundance(report, read_length=150)

TODO:
    - [ ] Add confidence score filtering
    - [ ] Implement batch processing for multiple samples
    - [ ] Add taxonomic summary table generation
"""

import logging
import subprocess
from pathlib import Path

logger = logging.getLogger(__name__)


def classify_reads(
    fastq: str,
    database: str,
    output_dir: str,
    threads: int = 8,
    confidence: float = 0.0,
) -> dict:
    """Classify metagenomic reads using Kraken2.

    Args:
        fastq: Path to input FASTQ file (single-end or interleaved).
        database: Path to Kraken2 database directory.
        output_dir: Directory for classification output and reports.
        threads: Number of threads for classification.
        confidence: Confidence threshold (0.0-1.0) for classification.

    Returns:
        Dictionary with paths to output file, report, and summary stats.

    Raises:
        FileNotFoundError: If input file or database doesn't exist.
        RuntimeError: If Kraken2 fails.
    """
    fastq_path = Path(fastq)
    db_path = Path(database)
    out_path = Path(output_dir)
    out_path.mkdir(parents=True, exist_ok=True)

    if not fastq_path.exists():
        raise FileNotFoundError(f"FASTQ not found: {fastq}")
    if not db_path.exists():
        raise FileNotFoundError(f"Database not found: {database}")

    sample_id = fastq_path.stem.replace(".fastq", "")
    output_file = out_path / f"{sample_id}.kraken2.output"
    report_file = out_path / f"{sample_id}.kraken2.report"

    cmd = [
        "kraken2",
        "--db", str(db_path),
        "--output", str(output_file),
        "--report", str(report_file),
        "--threads", str(threads),
        "--confidence", str(confidence),
        str(fastq_path),
    ]

    logger.info(f"Running Kraken2 on {sample_id}...")
    result = subprocess.run(cmd, capture_output=True, text=True)

    if result.returncode != 0:
        raise RuntimeError(f"Kraken2 failed: {result.stderr[:500]}")

    # Parse summary from stderr
    stats = _parse_kraken2_summary(result.stderr)

    logger.info(
        f"Kraken2 complete: {stats.get('classified_pct', 'N/A')} classified"
    )

    return {
        "output": str(output_file),
        "report": str(report_file),
        "stats": stats,
    }


def estimate_abundance(
    kraken_report: str,
    database: str,
    output_dir: str,
    read_length: int = 150,
    taxonomic_level: str = "S",
) -> dict:
    """Estimate species abundance using Bracken.

    Re-estimates abundance from Kraken2 report using Bayesian
    redistribution of reads assigned to higher taxonomic levels.

    Args:
        kraken_report: Path to Kraken2 report file.
        database: Path to Kraken2/Bracken database directory.
        output_dir: Directory for Bracken output.
        read_length: Read length for kmer length selection (50, 100, 150, 200, 250).
        taxonomic_level: Taxonomic level for abundance estimation.
            S=Species, G=Genus, F=Family, O=Order, C=Class, P=Phylum.

    Returns:
        Dictionary with path to abundance table and top species.
    """
    report_path = Path(kraken_report)
    out_path = Path(output_dir)
    out_path.mkdir(parents=True, exist_ok=True)

    sample_id = report_path.stem.replace(".kraken2.report", "")
    output_file = out_path / f"{sample_id}.bracken.{taxonomic_level}.tsv"

    cmd = [
        "bracken",
        "-d", str(database),
        "-i", str(report_path),
        "-o", str(output_file),
        "-r", str(read_length),
        "-l", taxonomic_level,
    ]

    logger.info(f"Running Bracken at level {taxonomic_level}...")
    result = subprocess.run(cmd, capture_output=True, text=True)

    if result.returncode != 0:
        raise RuntimeError(f"Bracken failed: {result.stderr[:500]}")

    return {
        "abundance_table": str(output_file),
        "taxonomic_level": taxonomic_level,
    }


def _parse_kraken2_summary(stderr_text: str) -> dict:
    """Parse Kraken2 classification summary from stderr output.

    Args:
        stderr_text: Raw stderr from Kraken2 run.

    Returns:
        Dictionary with total_reads, classified, unclassified, and percentages.
    """
    stats = {}
    for line in stderr_text.strip().split("\n"):
        line = line.strip()
        if "sequences classified" in line.lower():
            parts = line.split()
            for i, part in enumerate(parts):
                if part.startswith("("):
                    stats["classified_pct"] = part.strip("()")
        elif "sequences unclassified" in line.lower():
            parts = line.split()
            for i, part in enumerate(parts):
                if part.startswith("("):
                    stats["unclassified_pct"] = part.strip("()")
        elif "processed" in line.lower():
            parts = line.split()
            try:
                stats["total_reads"] = int(parts[0])
            except (ValueError, IndexError):
                pass

    return stats
