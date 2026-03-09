"""Support stage: secondary-structure features for shortlist tie-breaking.

This module is not one of the five core stations, but it gives us extra context
for The Winning Bunch when two candidates have similar enrichment behavior.
"""

import logging
import re
from dataclasses import dataclass

logger = logging.getLogger("aptamer_pipeline")

# We prefer ViennaRNA for stronger thermodynamic estimates.
# If it is not installed, we keep structure fields explicitly neutral so the
# rest of the pipeline does not mistake placeholders for real thermodynamics.
try:
    import RNA
    VIENNA_AVAILABLE = True
except ImportError:
    VIENNA_AVAILABLE = False
    logger.warning(
        "ViennaRNA not installed. Structure outputs will stay neutral. "
        "Tip: install with 'pip install ViennaRNA' for thermodynamic scoring."
    )


@dataclass
class StructureResult:
    """Secondary structure prediction result."""
    aptamer_id: str
    sequence: str
    structure: str  # dot-bracket notation
    mfe: float      # minimum free energy (kcal/mol)
    n_stems: int
    n_loops: int
    n_bulges: int
    has_g_quadruplex: bool
    motif_count: int

    def to_dict(self) -> dict:
        return {
            "aptamer_id": self.aptamer_id,
            "structure": self.structure,
            "mfe": round(self.mfe, 2),
            "n_stems": self.n_stems,
            "n_loops": self.n_loops,
            "n_bulges": self.n_bulges,
            "has_g_quadruplex": self.has_g_quadruplex,
            "motif_count": self.motif_count,
        }


def detect_g_quadruplex(sequence: str) -> bool:
    """Detect potential G-quadruplex forming sequences (G3+N1-7G3+N1-7G3+N1-7G3+)."""
    pattern = r"G{3,}[ACGT]{1,7}G{3,}[ACGT]{1,7}G{3,}[ACGT]{1,7}G{3,}"
    return bool(re.search(pattern, sequence.upper()))


def count_structural_motifs(structure: str) -> dict:
    """Count stems, loops, and bulges from dot-bracket notation."""
    n_stems = 0
    n_loops = 0
    n_bulges = 0

    i = 0
    while i < len(structure):
        if structure[i] == "(":
            # Consecutive paired bases are treated as one stem block.
            j = i
            while j < len(structure) and structure[j] == "(":
                j += 1
            n_stems += 1
            i = j
        elif structure[i] == ".":
            # Unpaired region; we classify it as loop or bulge using local context.
            j = i
            while j < len(structure) and structure[j] == ".":
                j += 1
            # If a dot-run sits between paired regions, we treat it as loop/bulge.
            if i > 0 and j < len(structure):
                if structure[i - 1] in "()" and structure[j] in "()":
                    if j - i <= 3:
                        n_bulges += 1
                    else:
                        n_loops += 1
                else:
                    n_loops += 1
            i = j
        else:
            i += 1

    return {"n_stems": n_stems, "n_loops": n_loops, "n_bulges": n_bulges}


def predict_structure(aptamer_id: str, sequence: str,
                      temperature: float = 37.0) -> StructureResult:
    """Predict secondary structure for a single aptamer sequence."""
    if VIENNA_AVAILABLE:
        md = RNA.md()
        md.temperature = temperature
        fc = RNA.fold_compound(sequence, md)
        structure, mfe = fc.mfe()
    has_gq = detect_g_quadruplex(sequence)
    if VIENNA_AVAILABLE:
        motifs = count_structural_motifs(structure)
        motif_count = motifs["n_stems"] + motifs["n_loops"] + motifs["n_bulges"]
        if has_gq:
            motif_count += 1
    else:
        structure = "NA"
        mfe = 0.0
        motifs = {"n_stems": 0, "n_loops": 0, "n_bulges": 0}
        motif_count = 0

    return StructureResult(
        aptamer_id=aptamer_id,
        sequence=sequence,
        structure=structure,
        mfe=mfe,
        n_stems=motifs["n_stems"],
        n_loops=motifs["n_loops"],
        n_bulges=motifs["n_bulges"],
        has_g_quadruplex=has_gq,
        motif_count=motif_count,
    )


def predict_structures(candidates: list, config: dict) -> list[StructureResult]:
    """Predict structures for all aptamer candidates.

    Args:
        candidates: List of AptamerCandidate objects.
        config: Full pipeline config dict.

    Returns:
        List of StructureResult objects.
    """
    temperature = config.get("structure", {}).get("temperature", 37.0)
    logger.info(f"Predicting structures for {len(candidates)} candidates "
                f"(T={temperature}°C, ViennaRNA={'yes' if VIENNA_AVAILABLE else 'no'})")

    results = []
    for candidate in candidates:
        result = predict_structure(candidate.id, candidate.sequence, temperature)
        results.append(result)

    logger.info(f"Structure prediction complete. "
                f"MFE range: {min(r.mfe for r in results):.1f} to "
                f"{max(r.mfe for r in results):.1f} kcal/mol")

    return results
