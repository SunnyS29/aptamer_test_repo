# Empirical HT-SELEX Winner Pipeline

This project helps us move from raw HT-SELEX count tables to a ranked shortlist of aptamer candidates.
The key idea is simple: we reward sequences that truly gain ground across rounds, not sequences that only look big in one file.
Everything in this README is written so a new teammate can explain both the workflow and the failure points with confidence.

## What This Pipeline Is (and Is Not)

- **It is** an empirical analysis pipeline that reads real sequence counts from SELEX rounds.
- **It is not** a random-sequence generator.
- **It is designed** to make our decisions traceable: we can explain why each sequence was kept, filtered, or ranked.

## The 5 Stations

### 1. Scanner (Data Ingestion)
- We read a counts table and detect sequence + round columns.
- We support both wide format (`round_1`, `round_2`, ...) and long format (`round`, `count`).
- We merge duplicate sequence rows and keep clean DNA-style sequence strings.

### 2. Race Begins (CPM Normalization)
- We compute each round's total reads.
- We convert raw counts to CPM (`count / round_total * 1,000,000`) so rounds with different depths are comparable.
- This is what lets us interpret "growth" as a biological trend instead of a sequencing-depth artifact.

### 3. Winning Bunch (Enrichment Scoring)
- For each sequence, we look at how CPM changes across rounds.
- We use two signals together: first-to-last log2 enrichment and overall trend slope.
- Sequences with strong and steady rise get higher enrichment scores.

### 4. Security Check (Target Verification)
- We fetch target information from PDB/UniProt (or read FASTA).
- If target retrieval fails, we hard-stop the run.
- We do this to protect us from silent junk outputs based on missing target context.

### 5. Final Cut (Filtering + Ranking)
- We remove weak candidates using enrichment and structure thresholds.
- We compute a composite score where enrichment is the primary driver, and structure/diversity help break ties.
- The output shortlist is saved to CSV/JSON for downstream review.

## Quick Start

### Install
```bash
git clone https://github.com/SunnyS29/aptamer_test_repo.git
cd aptamer_test_repo
pip install -r requirements.txt
```

Optional for stronger structure modeling:
```bash
pip install ViennaRNA
```

### Run
```bash
python -m src.pipeline --config config/pipeline_config.yaml
```

Run one station for debugging:
```bash
python -m src.pipeline --config config/pipeline_config.yaml --stage library
```

## Input Format We Expect

### Wide format
```csv
sequence,round_1,round_2,round_3
ACGT...,10,25,80
TGCA...,8,12,9
```

### Long format
```csv
sequence,round,count
ACGT...,round_1,10
ACGT...,round_2,25
```

## Configuration Cheat Sheet

Edit `config/pipeline_config.yaml`:

- `target`: where target info comes from (`pdb_id`, `fasta`, `smiles`, `uniprot`)
- `selex.counts_file`: the counts table path
- `library`: sequence QC filters (length, GC, homopolymer, min total count)
- `scoring`: pseudocount + growth weights
- `filtering`: shortlist strictness (`top_n`, `min_log2_enrichment`, structure limits)
- `output`: file format + output directory

## Friendly Troubleshooting (By Station)

### Scanner
- If the pipeline says **"Could not find a sequence column"**, don't panic.
- It usually means column headers are inconsistent.
- Tip: rename the sequence column to `sequence` and rerun.

### Race Begins
- If you see **"One or more rounds have zero total reads"**, normalization cannot proceed safely.
- This means at least one round is effectively empty.
- Tip: check that all round columns are mapped correctly and not accidentally blank.

### Winning Bunch
- If scores look flat (many near 0.5), your trajectories may be too similar or too sparse.
- That is not a code crash, but it is a signal to inspect round quality and `min_total_count`.
- Tip: inspect `log2_enrichment` and `trend_slope` in exported results.

### Security Check
- If you see **"Failed to fetch PDB target"** or **"No sequence found"**, the run is correctly blocking unsafe analysis.
- Tip: check network access, ID spelling, or switch to a local FASTA target file.
- We prefer stopping early here because silent fallback would corrupt later decisions.

### Final Cut
- If you get **"No candidates passed filters"**, the filters are likely too strict for the current dataset.
- Tip: loosen `min_log2_enrichment`, `mfe_threshold`, or `min_complexity` gradually and rerun.
- Keep a record of threshold changes so we can justify shortlist criteria later.

## Project Structure

```text
src/
├── pipeline.py            # Orchestrates the full run and stage-by-stage execution
├── sequence_generator.py  # Scanner + Race Begins logic
├── binding_scorer.py      # Winning Bunch scoring
├── target_analyzer.py     # Security Check target validation
├── filter_rank.py         # Final Cut filtering and ranking
├── structure_predictor.py # Secondary structure features
└── utils.py               # Shared helpers
```

## Testing

```bash
python -m pytest tests/ -v
```
