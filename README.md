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

## How The Race Works

This is the plain-language version of what the pipeline is doing under the hood.
Think of every sequence as a runner in a stadium. We are not just asking who looked good once. We are asking who keeps moving toward the front as the race gets harder.

### 1. The Trimmer
- Job: if raw FASTQ reads still contain constant primer regions, we cut those away and keep only the variable insert.
- Why it matters: otherwise we risk ranking sequencing constructs instead of real aptamer candidates.
- Math behind it: pattern matching, not a statistical model. We look for a left anchor and a right anchor and keep the sequence between them. We can also check the reverse complement if reads come in the opposite orientation.

### 2. The Counter
- Job: count how many times each sequence appears in each round.
- Why it matters: this gives us the raw race table.
- Math behind it: empirical frequency counting. No prediction yet, just observed abundance.

### 3. The Equalizer
- Job: convert raw counts into CPM, or counts per million.
- Formula: `CPM = count / total_round_reads * 1,000,000`
- Why it matters: rounds often have different sequencing depth, so raw counts alone are not fair.
- Simple meaning: CPM tells us how much of the pool a sequence owns in that round.

### 4. The Growth Judge
- Job: score how strongly and how steadily a sequence grows across rounds.
- Exact methods used:
- `log2` enrichment: measures doubling-like growth from first round to last round.
- Least-squares slope: fits a straight trend line across rounds to measure overall upward movement.
- Monotonicity: checks how often a sequence avoids going backwards from one round to the next.
- RMSE: measures how far the real trajectory wiggles away from the fitted trend line.
- Min-max scaling: rescales the growth features to `0-1` so they can be combined fairly.
- Weighted sum: fold change `0.60`, slope `0.20`, pace consistency `0.20`.
- Simple meaning: we reward sequences that climb hard and climb smoothly.

### 5. The Shape Check
- Job: add secondary-structure context as a supporting signal.
- Exact methods used:
- ViennaRNA minimum free energy, if installed.
- Dot-bracket motif counting for stems, loops, and bulges.
- G-quadruplex pattern detection with a sequence regex.
- Important note: if ViennaRNA is not installed, structure now stays neutral instead of pretending the fallback is real biology.

### 6. The Winning Bunch
- Job: build the final shortlist.
- Exact methods used:
- Stability score: min-max normalized MFE.
- Diversity score: k-mer rarity across the candidate pool.
- Weighted sum: enrichment growth `0.70`, structural stability `0.20`, sequence diversity `0.10`.
- Simple meaning: growth is the main decision-maker, while structure and diversity help break ties.

### 7. The Finish-Line Referee
- Job: decide whether the experiment looks mature enough to stop.
- Exact methods used:
- Top-10 overlap percentage between the last two rounds.
- Jaccard index for set similarity of the two leaderboards.
- Coverage percentages for the top 1, top 10, and top 100 sequences.
- Mean acceleration of the top 3 trajectories.
- Mean, median, and coefficient of variation of the pace score.
- A composite data quality score and a rule-based recommendation.

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

### Convert Raw Sequencing Files (Simple 3-Step Guide)
Use this if you have one file per round (`.fastq`, `.fastq.gz`, `.fasta`, or `.fasta.gz`).

1. Put your round files in one folder.
2. Run this command (copy/paste, then replace file names):

```bash
python -m src.fasta_round_counter \
  data/PRJDB9110/DRR201861.fastq.gz \
  data/PRJDB9110/DRR201862.fastq.gz \
  data/PRJDB9110/DRR201863.fastq.gz \
  --output data/PRJDB9110/prjdb9110_round_counts.csv \
  --round-labels round_0 round_1 round_2 \
  --summary data/PRJDB9110/prjdb9110_round_summary.tsv
```

3. Tell the pipeline to use the new counts file:

```yaml
selex:
  counts_file: "data/PRJDB9110/prjdb9110_round_counts.csv"
```

What this does:
- Reads each round file.
- Counts how many times each sequence appears.
- Builds one table the pipeline can compare across rounds.

Helpful tip:
- If you do not pass `--round-labels`, we try to guess rounds from file names like `round_3`, `r3`, or `rnd3`.
- If round numbers are not in file names, files are sorted alphabetically.
- If your FASTQ reads still include constant primer regions, add `--left-anchor` and `--right-anchor` so we count only the variable insert.

Example with anchor extraction:
```bash
python -m src.fasta_round_counter \
  data/run_round1.fastq.gz \
  data/run_round2.fastq.gz \
  --left-anchor AGACGCAACTGAATGAA \
  --right-anchor CCGTAACTAGTCGCGTCAC \
  --output data/run_counts.csv
```

## Supported Input Formats

Use any one of these options.

### Option 1 (Easiest): One table with one row per sequence
Your file can be `.csv` or `.tsv`.

```csv
sequence,round_1,round_2,round_3
ACGT...,10,25,80
TGCA...,8,12,9
```

What matters:
- First column is sequence text (`sequence`).
- Each round has its own count column (`round_1`, `round_2`, etc.).

### Option 2: Long table (one row per sequence per round)
Use this if your export is already in long format.

```csv
sequence,round,count
ACGT...,round_1,10
ACGT...,round_2,25
```

What matters:
- Must include all three columns: `sequence`, `round`, `count`.
- Counts must be whole numbers (no decimals).

### Option 3: Raw sequencing files (FASTQ/FASTA)
Use this when you start from raw files from the sequencer.

- Supported: `.fastq`, `.fastq.gz`, `.fasta`, `.fasta.gz`
- One file should represent one round.
- Convert first with `python -m src.fasta_round_counter ...`
- Then run the pipeline on the new counts table.

Common input mistakes:
- Missing `sequence` column name.
- A round column with all zeros.
- Mixed files from different experiments in one run.

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
- If ViennaRNA is not installed, structure stays neutral instead of steering the ranking.

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
