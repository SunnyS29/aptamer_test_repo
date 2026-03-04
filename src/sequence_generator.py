"""Stage 2: SELEX-inspired ssDNA aptamer library generation."""

import random
import logging
from dataclasses import dataclass
from src.utils import gc_content, has_homopolymer

logger = logging.getLogger("aptamer_pipeline")

NUCLEOTIDES = "ACGT"


@dataclass
class AptamerCandidate:
    """Represents a generated aptamer candidate."""
    id: str
    sequence: str
    length: int
    gc: float

    def to_dict(self) -> dict:
        return {
            "id": self.id,
            "sequence": self.sequence,
            "length": self.length,
            "gc_content": round(self.gc, 4),
        }


def generate_random_sequence(length: int, rng: random.Random) -> str:
    """Generate a random DNA sequence of given length."""
    return "".join(rng.choice(NUCLEOTIDES) for _ in range(length))


def validate_sequence(seq: str, gc_min: float, gc_max: float,
                      max_homopolymer: int) -> bool:
    """Check if a sequence passes quality filters."""
    gc = gc_content(seq)
    if gc < gc_min or gc > gc_max:
        return False
    if has_homopolymer(seq, max_homopolymer):
        return False
    return True


def generate_library(config: dict) -> list[AptamerCandidate]:
    """Generate a filtered library of ssDNA aptamer candidates.

    Args:
        config: Dictionary with keys from the 'library' section of pipeline config.

    Returns:
        List of AptamerCandidate objects that pass all quality filters.
    """
    lib_config = config["library"]
    size = lib_config["size"]
    length_min = lib_config["length_min"]
    length_max = lib_config["length_max"]
    gc_min = lib_config["gc_min"]
    gc_max = lib_config["gc_max"]
    max_homo = lib_config["max_homopolymer"]
    seed = lib_config.get("seed", 42)

    rng = random.Random(seed)
    candidates = []
    attempts = 0
    max_attempts = size * 20  # safety limit

    logger.info(f"Generating aptamer library: target size={size}, "
                f"length={length_min}-{length_max}nt, "
                f"GC={gc_min}-{gc_max}")

    while len(candidates) < size and attempts < max_attempts:
        attempts += 1
        length = rng.randint(length_min, length_max)
        seq = generate_random_sequence(length, rng)

        if validate_sequence(seq, gc_min, gc_max, max_homo):
            apt_id = f"APT_{len(candidates) + 1:06d}"
            gc = gc_content(seq)
            candidates.append(AptamerCandidate(
                id=apt_id, sequence=seq, length=length, gc=gc
            ))

    logger.info(f"Generated {len(candidates)} candidates from {attempts} attempts "
                f"({len(candidates)/attempts*100:.1f}% pass rate)")

    return candidates
