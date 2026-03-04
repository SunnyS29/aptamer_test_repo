"""Shared utility functions for the aptamer pipeline."""

import logging
import os
import yaml
from pathlib import Path


def setup_logging(verbose: bool = True) -> logging.Logger:
    """Configure pipeline logging."""
    logger = logging.getLogger("aptamer_pipeline")
    logger.setLevel(logging.DEBUG if verbose else logging.INFO)

    if not logger.handlers:
        handler = logging.StreamHandler()
        handler.setFormatter(
            logging.Formatter("%(asctime)s [%(levelname)s] %(message)s",
                              datefmt="%Y-%m-%d %H:%M:%S")
        )
        logger.addHandler(handler)

    return logger


def load_config(config_path: str) -> dict:
    """Load pipeline configuration from YAML file."""
    path = Path(config_path)
    if not path.exists():
        raise FileNotFoundError(f"Config file not found: {config_path}")

    with open(path) as f:
        return yaml.safe_load(f)


def ensure_output_dir(output_dir: str) -> Path:
    """Create output directory if it doesn't exist."""
    path = Path(output_dir)
    path.mkdir(parents=True, exist_ok=True)
    return path


def gc_content(sequence: str) -> float:
    """Calculate GC content of a nucleotide sequence."""
    seq = sequence.upper()
    gc = sum(1 for nt in seq if nt in "GC")
    return gc / len(seq) if seq else 0.0


def has_homopolymer(sequence: str, max_run: int = 4) -> bool:
    """Check if sequence contains a homopolymer run exceeding max_run."""
    count = 1
    for i in range(1, len(sequence)):
        if sequence[i] == sequence[i - 1]:
            count += 1
            if count > max_run:
                return True
        else:
            count = 1
    return False


def read_fasta(filepath: str) -> list[tuple[str, str]]:
    """Read sequences from a FASTA file. Returns list of (header, sequence) tuples."""
    sequences = []
    header = None
    seq_parts = []

    with open(filepath) as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if header is not None:
                    sequences.append((header, "".join(seq_parts)))
                header = line[1:].strip()
                seq_parts = []
            elif line:
                seq_parts.append(line.upper())

    if header is not None:
        sequences.append((header, "".join(seq_parts)))

    return sequences


def write_fasta(sequences: list[tuple[str, str]], filepath: str) -> None:
    """Write sequences to a FASTA file."""
    with open(filepath, "w") as f:
        for header, seq in sequences:
            f.write(f">{header}\n")
            for i in range(0, len(seq), 80):
                f.write(seq[i:i + 80] + "\n")
