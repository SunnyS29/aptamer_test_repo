# Aptamer Target Identification Pipeline

A computational pipeline for accelerating DNA aptamer target identification, inspired by electrochemical aptamer-based (EAB) biosensor applications such as wearable drug monitoring patches.

## Overview

This pipeline generates, evaluates, and ranks ssDNA aptamer candidates for a given molecular target using:

1. **Target Analysis** — Parse protein/small molecule targets, extract biochemical features
2. **Library Generation** — SELEX-inspired random ssDNA aptamer library with quality filters
3. **Structure Prediction** — Secondary structure prediction via ViennaRNA (with fallback estimator)
4. **Binding Scoring** — ML-based binding affinity prediction (Random Forest / Gradient Boosting)
5. **Filtering & Ranking** — Multi-criteria composite scoring and candidate selection
6. **Results Export** — CSV/JSON output with summary visualizations

## Installation

```bash
git clone https://github.com/SunnyS29/aptamer_test_repo.git
cd aptamer_test_repo
pip install -r requirements.txt
```

Optional (for accurate structure prediction):
```bash
pip install ViennaRNA
```

## Usage

Run the full pipeline:
```bash
python -m src.pipeline --config config/pipeline_config.yaml
```

Run a specific stage:
```bash
python -m src.pipeline --config config/pipeline_config.yaml --stage library
```

Override top N candidates:
```bash
python -m src.pipeline --config config/pipeline_config.yaml --top-n 100
```

## Configuration

Edit `config/pipeline_config.yaml` to customize:

- **target**: Input type (PDB ID, FASTA, SMILES, UniProt), target name
- **library**: Size, length range, GC content bounds, homopolymer limits
- **structure**: Folding temperature, MFE threshold
- **scoring**: ML model type, number of estimators, scoring weights
- **filtering**: Top N candidates, minimum structural complexity
- **output**: Format (CSV/JSON), plot generation

## Project Structure

```
src/
├── pipeline.py            # CLI orchestrator
├── sequence_generator.py  # Aptamer library generation
├── target_analyzer.py     # Target molecule analysis
├── structure_predictor.py # Secondary structure prediction
├── binding_scorer.py      # ML binding affinity scoring
├── filter_rank.py         # Filtering and ranking
└── utils.py               # Shared utilities
```

## Running Tests

```bash
pip install pytest
python -m pytest tests/ -v
```

## Notes

- The ML scoring model uses synthetic training data as a baseline. Replace with experimental SELEX data for production use.
- ViennaRNA provides accurate thermodynamic structure prediction. Without it, a simplified estimator is used.
- G-quadruplex detection is included as these motifs are common in functional aptamers.
