"""Stage 1: Target molecule input and feature extraction."""

import logging
import requests
from dataclasses import dataclass, field

logger = logging.getLogger("aptamer_pipeline")

# Amino acid hydrophobicity scale (Kyte-Doolittle)
HYDROPHOBICITY = {
    "I": 4.5, "V": 4.2, "L": 3.8, "F": 2.8, "C": 2.5,
    "M": 1.9, "A": 1.8, "G": -0.4, "T": -0.7, "S": -0.8,
    "W": -0.9, "Y": -1.3, "P": -1.6, "H": -3.2, "D": -3.5,
    "E": -3.5, "N": -3.5, "Q": -3.5, "K": -3.9, "R": -4.5,
}

# Amino acid charge at pH 7.4
CHARGE = {
    "D": -1, "E": -1, "K": 1, "R": 1, "H": 0.1,
}


@dataclass
class TargetFeatures:
    """Extracted features of the target molecule."""
    name: str
    input_type: str
    sequence: str = ""
    molecular_weight: float = 0.0
    avg_hydrophobicity: float = 0.0
    net_charge: float = 0.0
    length: int = 0
    aromatic_fraction: float = 0.0
    polar_fraction: float = 0.0
    metadata: dict = field(default_factory=dict)

    def to_dict(self) -> dict:
        return {
            "name": self.name,
            "input_type": self.input_type,
            "length": self.length,
            "molecular_weight": round(self.molecular_weight, 2),
            "avg_hydrophobicity": round(self.avg_hydrophobicity, 4),
            "net_charge": round(self.net_charge, 2),
            "aromatic_fraction": round(self.aromatic_fraction, 4),
            "polar_fraction": round(self.polar_fraction, 4),
        }

    def to_feature_vector(self) -> list[float]:
        """Return numeric features for ML scoring."""
        return [
            self.molecular_weight,
            self.avg_hydrophobicity,
            self.net_charge,
            self.length,
            self.aromatic_fraction,
            self.polar_fraction,
        ]


def fetch_pdb_sequence(pdb_id: str) -> str:
    """Fetch protein sequence from RCSB PDB."""
    url = f"https://data.rcsb.org/rest/v1/core/polymer_entity/{pdb_id}/1"
    logger.info(f"Fetching PDB entry: {pdb_id}")

    try:
        resp = requests.get(url, timeout=15)
        resp.raise_for_status()
        data = resp.json()
        sequence = data.get("entity_poly", {}).get("pdbx_seq_one_letter_code_can", "")
        if sequence:
            logger.info(f"Retrieved sequence: {len(sequence)} residues")
            return sequence
        raise ValueError(f"No sequence found for PDB ID {pdb_id}")
    except requests.RequestException as e:
        logger.warning(f"PDB fetch failed: {e}. Using offline mode.")
        return ""


def fetch_uniprot_sequence(uniprot_id: str) -> str:
    """Fetch protein sequence from UniProt."""
    url = f"https://rest.uniprot.org/uniprotkb/{uniprot_id}.fasta"
    logger.info(f"Fetching UniProt entry: {uniprot_id}")

    try:
        resp = requests.get(url, timeout=15)
        resp.raise_for_status()
        lines = resp.text.strip().split("\n")
        sequence = "".join(line for line in lines if not line.startswith(">"))
        logger.info(f"Retrieved sequence: {len(sequence)} residues")
        return sequence
    except requests.RequestException as e:
        logger.warning(f"UniProt fetch failed: {e}. Using offline mode.")
        return ""


def compute_features(sequence: str) -> dict:
    """Compute biochemical features from an amino acid sequence."""
    if not sequence:
        return {}

    seq = sequence.upper()
    length = len(seq)

    # Average molecular weight per residue ~110 Da
    mw = length * 110.0

    # Hydrophobicity
    hydro_scores = [HYDROPHOBICITY.get(aa, 0.0) for aa in seq]
    avg_hydro = sum(hydro_scores) / length if length else 0.0

    # Net charge
    net_charge = sum(CHARGE.get(aa, 0.0) for aa in seq)

    # Aromatic fraction (F, W, Y)
    aromatic = sum(1 for aa in seq if aa in "FWY") / length

    # Polar fraction (S, T, N, Q, D, E, K, R, H)
    polar = sum(1 for aa in seq if aa in "STNQDEKRH") / length

    return {
        "molecular_weight": mw,
        "avg_hydrophobicity": avg_hydro,
        "net_charge": net_charge,
        "aromatic_fraction": aromatic,
        "polar_fraction": polar,
    }


def analyze_target(config: dict) -> TargetFeatures:
    """Analyze the target molecule and extract features.

    Args:
        config: Full pipeline config dict.

    Returns:
        TargetFeatures with extracted biochemical properties.
    """
    target_config = config["target"]
    input_type = target_config["input_type"]
    input_value = target_config["input_value"]
    name = target_config.get("name", input_value)

    logger.info(f"Analyzing target: {name} (type={input_type})")

    sequence = ""
    metadata = {}

    if input_type == "pdb_id":
        sequence = fetch_pdb_sequence(input_value)
        metadata["pdb_id"] = input_value
    elif input_type == "fasta":
        from src.utils import read_fasta
        seqs = read_fasta(input_value)
        if seqs:
            sequence = seqs[0][1]
            metadata["fasta_header"] = seqs[0][0]
    elif input_type == "smiles":
        # For small molecules, store SMILES directly
        metadata["smiles"] = input_value
        logger.info("Small molecule target (SMILES) — limited feature extraction")
        return TargetFeatures(
            name=name, input_type=input_type, sequence=input_value,
            metadata=metadata
        )
    elif input_type == "uniprot":
        sequence = fetch_uniprot_sequence(input_value)
        metadata["uniprot_id"] = input_value

    features = compute_features(sequence)

    return TargetFeatures(
        name=name,
        input_type=input_type,
        sequence=sequence,
        length=len(sequence),
        molecular_weight=features.get("molecular_weight", 0.0),
        avg_hydrophobicity=features.get("avg_hydrophobicity", 0.0),
        net_charge=features.get("net_charge", 0.0),
        aromatic_fraction=features.get("aromatic_fraction", 0.0),
        polar_fraction=features.get("polar_fraction", 0.0),
        metadata=metadata,
    )
