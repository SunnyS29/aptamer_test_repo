# Aptamer Playground Race:
## An Empirical HT-SELEX Enrichment & Analysis Pipeline

This project helps us move from raw HT-SELEX count tables to a ranked shortlist of aptamer candidates.
The key idea is simple: we reward sequences that truly gain ground across rounds, not sequences that only look big in one file.
Everything in this README is written so that anyone can interact with the tool and be aware of the processes behind generating the output.

## What This Pipeline Is (and Is Not)

- **It is** an empirical analysis pipeline that reads real sequence counts from SELEX rounds.
- **It is not** a random-sequence generator.
- **It is designed** to make our decisions traceable: we can explain why each sequence was kept, filtered, or ranked.

## The 5 Stations

### 1. The Scanner (Data Ingestion)
- We read a counts table and detect sequence + round columns.
- We support both wide format (`round_1`, `round_2`, ...) and long format (`round`, `count`).
- We merge duplicate sequence rows and keep clean DNA-style sequence strings.

### 2. The Starting Line (CPM Normalization)
- We compute each round's total reads.
- We convert raw counts to CPM (`count / round_total * 1,000,000`) so rounds with different depths are comparable.
- This is what lets us interpret "growth" as a biological trend instead of a sequencing-depth artifact.

### 3. The Race Begins (Enrichment Scoring)
- For each sequence, we look at how CPM changes across rounds.
- We use three signals together: first-to-last log2 enrichment, overall trend slope, and pace consistency.
- Pace consistency rewards sequences that climb smoothly instead of jumping late in one spike.
- Sequences with strong and steady rise get higher enrichment scores.

### 4. Security Check (Target Verification)
- We fetch target information from PDB/UniProt (or read FASTA).
- If target retrieval fails, we hard-stop the run.
- We do this to protect us from silent junk outputs based on missing target context.

### 5. The Winning Bunch (Filtering + Ranking)
- We remove weak candidates using enrichment and structure thresholds.
- We compute a composite score where enrichment is the primary driver, and structure/diversity help break ties.
- Diversity is computed with a pooled k-mer rarity score (near-linear runtime), not all-vs-all pairwise distances.
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

### Convert FASTQ/FASTA Rounds to Pipeline Input
If your lab gives you one sequencing file per round (`.fastq`, `.fastq.gz`, `.fasta`, or `.fasta.gz`), run this first:

```bash
python -m src.fasta_round_counter \
  data/PRJDB9110/DRR201861.fastq.gz \
  data/PRJDB9110/DRR201862.fastq.gz \
  data/PRJDB9110/DRR201863.fastq.gz \
  --output data/PRJDB9110/prjdb9110_round_counts.csv \
  --round-labels round_0 round_1 round_2 \
  --summary data/PRJDB9110/prjdb9110_round_summary.tsv
```

Then point your config to the new counts table:

```yaml
selex:
  counts_file: "data/PRJDB9110/prjdb9110_round_counts.csv"
```

Tip: if you skip `--round-labels`, we infer round names from filenames like `round_3`, `r3`, or `rnd3`.
If names do not include round numbers, files are sorted alphabetically.

## Supported Input Formats

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

### Raw Round Files (Pre-Step)
- Supported: `.fastq`, `.fastq.gz`, `.fasta`, `.fasta.gz`
- Each file should represent one SELEX round.
- Use `python -m src.fasta_round_counter ...` to build the wide table above.

## Configuration Cheat Sheet

Edit `config/pipeline_config.yaml`:

- `target`: where target info comes from (`pdb_id`, `fasta`, `smiles`, `uniprot`)
- `selex.counts_file`: the counts table path
- `library`: sequence QC filters (length, GC, homopolymer, min total count)
- `scoring`: pseudocount + growth weights
- `scoring.vectorized_metrics`: set `true` to speed up pace/slope calculations with NumPy on large libraries (default `false`)
- `scoring.diversity_kmer_size`: k-mer size used for diversity rarity scoring (default `3`)
- `filtering`: shortlist strictness (`top_n`, `min_log2_enrichment`, structure limits)
- `output`: file format + output directory

## Friendly Troubleshooting (By Section)

### The Scanner
- If the pipeline says **"Could not find a sequence column"**, don't panic.
- It usually means column headers are inconsistent.
- Tip: rename the sequence column to `sequence` and rerun.
- If conversion fails with **"Unsupported input format"**, check file suffixes.
- Tip: rename files to `.fastq(.gz)` or `.fasta(.gz)` and rerun the converter.

### The Starting Line
- If you see **"One or more rounds have zero total reads"**, normalization cannot proceed safely.
- This means at least one round is effectively empty.
- Tip: check that all round columns are mapped correctly and not accidentally blank.

### The Race Begins
- If scores look flat (many near 0.5), your trajectories may be too similar or too sparse.
- That is not a code crash, but it is a signal to inspect round quality and `min_total_count`.
- Tip: inspect `log2_enrichment`, `trend_slope`, and `pace_consistency` in exported results.
- If this stage is slow with very large candidate sets, try `scoring.vectorized_metrics: true`.

### Security Check
- If you see **"Failed to fetch PDB target"** or **"No sequence found"**, the run is correctly blocking unsafe analysis.
- Tip: check network access, ID spelling, or switch to a local FASTA target file.
- We prefer stopping early here because silent fallback would corrupt later decisions.

### The Winning Bunch
- If you get **"No candidates passed filters"**, the filters are likely too strict for the current dataset.
- Tip: loosen `min_log2_enrichment`, `mfe_threshold`, or `min_complexity` gradually and rerun.
- Keep a record of threshold changes so we can justify shortlist criteria later.
- If you see **"scoring.diversity_kmer_size must be >= 1"**, set `scoring.diversity_kmer_size` to `3` and rerun.
- If ranking still feels slow, raise `library.min_total_count` to reduce the candidate pool before Station 5.

## Project Structure

```text
src/
├── pipeline.py            # Orchestrates the full run and stage-by-stage execution
├── sequence_generator.py  # The Scanner + The Starting Line logic
├── binding_scorer.py      # The Race Begins scoring
├── target_analyzer.py     # Security Check target validation
├── filter_rank.py         # The Winning Bunch filtering and ranking
├── structure_predictor.py # Secondary structure features
└── utils.py               # Shared helpers
```

## Testing

```bash
python -m pytest tests/ -v
```

## Stopping-Point Diagnostic

Use this when you want a quick health check on whether SELEX rounds are converging:

```bash
python -m src.stopping_diagnostic --config config/pipeline_config.yaml
```

It reports four universal markers:
- Leaderboard stability between the last two rounds
- Top-candidate slope trajectory (acceleration/deceleration)
- Pool dominance coverage (top 1 / 10 / 100)
- Pace consistency across top-ranked candidates

Recommendation output:
- `A`: Sequence more rounds
- `B`: Stop and validate
- `C`: Potential over-selection, review earlier rounds
