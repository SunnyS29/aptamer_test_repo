# Aptamer Pipeline
## An Empirical HT-SELEX Enrichment & Analysis Pipeline

This pipeline moves from raw HT-SELEX count tables to a ranked shortlist of aptamer candidates.
The key idea is simple: sequences are rewarded for truly gaining ground across rounds, not for looking big in one snapshot.
Everything in this README is written so that anyone can interact with the tool and understand the processes behind the output.

## What This Pipeline Is (and Is Not)

- **It is** an empirical analysis pipeline that reads real sequence counts from SELEX rounds.
- **It is not** a random-sequence generator.
- **It is designed** to make decisions traceable: every sequence that is kept, filtered, or ranked can be explained.
- **It is not** a reliable aptamer structure predictor. Rough structure annotations can be attached, but they are not used to rank winners.
- **It does not** replace wet-lab validation. The shortlist is meant to narrow the field, not prove binding on its own.

## The 5 Stations

### 1. The Scanner (Data Ingestion)
- Reads a counts table and detects sequence + round columns.
- Supports both wide format (`round_1`, `round_2`, ...) and long format (`round`, `count`).
- Merges duplicate sequence rows and keeps clean DNA-style sequence strings.

### 2. The Starting Line (CPM Normalisation)
- Computes each round's total reads.
- Converts raw counts to CPM (`count / round_total * 1,000,000`) so rounds with different depths are comparable.
- This is what allows "growth" to be interpreted as a biological trend rather than a sequencing-depth artifact.

### 3. The Race Begins (Enrichment Scoring)
- For each sequence, looks at how CPM changes across rounds.
- Uses three signals together: first-to-last log2 enrichment, overall trend slope, and a light terminal guardrail.
- Log2 enrichment and slope do the main ranking work.
- The terminal guardrail catches sequences that fade in the last round, without penalising step-function winners that take off late.

### 4. Security Check (Target Verification)
- Fetches target information from PDB/UniProt (or reads FASTA).
- If target retrieval fails, the run hard-stops.
- This protects against silent junk outputs based on missing target context.

### 5. The Winning Bunch (Filtering + Ranking)
- Removes weak candidates using enrichment thresholds.
- Computes a composite score where enrichment is the primary driver, and diversity is a light tie-breaker.
- Diversity is computed with a pooled k-mer rarity score (near-linear runtime), not all-vs-all pairwise distances.
- The output shortlist is saved to CSV/JSON for downstream review.

## How The Race Works

This is the plain-language version of what the pipeline is doing under the hood.
Think of every sequence as a runner in a stadium. The question is not who looked good once — it is who keeps moving toward the front as the race gets harder.

### 1. The Trimmer
- Job: if raw FASTQ reads still contain constant primer regions, cut those away and keep only the variable insert.
- Why it matters: otherwise the risk is ranking sequencing constructs instead of real aptamer candidates.
- Math behind it: pattern matching, not a statistical model. Looks for a left anchor and a right anchor and keeps the sequence between them. Can also check the reverse complement if reads come in the opposite orientation.

### 2. The Counter
- Job: count how many times each sequence appears in each round.
- Why it matters: this gives the raw race table.
- Math behind it: empirical frequency counting. No prediction yet, just observed abundance.

### 3. The Equalizer
- Job: convert raw counts into CPM, or counts per million.
- Formula: `CPM = count / total_round_reads * 1,000,000`
- Why it matters: rounds often have different sequencing depth, so raw counts alone are not fair.
- Simple meaning: CPM tells you how much of the pool a sequence owns in that round.

### 4. The Growth Judge
- Job: score how strongly a sequence takes over the pool, while checking that it does not fade at the finish.
- Exact methods used:
- `log2` enrichment: measures doubling-like growth from first round to last round.
- Least-squares slope: fits a straight trend line across rounds to measure overall upward movement.
- Monotonicity: checks how often a sequence avoids going backwards from one round to the next.
- RMSE: measures how far the real trajectory wiggles away from the fitted trend line.
- Terminal guardrail: checks whether the last round goes down versus the round before it.
- Min-max scaling: rescales the growth features to `0-1` so they can be combined fairly.
- Weighted sum: fold change `0.80`, slope `0.15`, terminal guardrail `0.05`.
- Simple meaning: sequences that take over the pool are rewarded, and candidates that fade at the end are lightly penalised.

### 5. The Shape Check
- Job: attach optional structure annotations for review.
- Exact methods used:
- ViennaRNA minimum free energy, if installed.
- Dot-bracket motif counting for stems, loops, and bulges.
- G-quadruplex pattern detection with a sequence regex.
- Important note: these fields are annotations only. They are not used to rank candidates because aptamer structure prediction is not reliable enough to drive winner selection.

### 6. The Winning Bunch
- Job: build the final shortlist.
- Exact methods used:
- Diversity score: k-mer rarity across the candidate pool.
- Weighted sum: enrichment growth `0.90`, sequence diversity `0.10`.
- Simple meaning: growth is the main decision-maker, while diversity helps break ties.

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

Optional for rough structure annotation only:
```bash
pip install ViennaRNA
```

### Run
```bash
python -m src.pipeline --config config/pipeline_config.yaml
```

If you want the tool to prompt for files instead of typing paths by hand:
```bash
python -m src.interactive_launcher
```

The launcher will:
- ask whether you already have a counts table or raw round files
- let you choose the target input
- let you choose an output folder
- save an `interactive_run_config.yaml` file so the run is still reproducible

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
- If you do not pass `--round-labels`, round numbers are guessed from file names like `round_3`, `r3`, or `rnd3`.
- If round numbers are not in file names, files are sorted alphabetically.
- If your FASTQ reads still include constant primer regions, add `--left-anchor` and `--right-anchor` to count only the variable insert.

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
Use this when starting from raw files from the sequencer.

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
- `scoring.growth_weights.terminal_guardrail`: optional name for the anti-fader weight; older configs may still use `pace_consistency`
- `scoring.diversity_kmer_size`: k-mer size used for diversity rarity scoring (default `3`)
- `filtering`: shortlist strictness (`top_n`, `min_log2_enrichment`)
- `output`: file format + output directory

## Friendly Troubleshooting (By Section)

### The Scanner
- If the pipeline says **"Could not find a sequence column"**, it usually means column headers are inconsistent.
- Tip: rename the sequence column to `sequence` and rerun.
- If conversion fails with **"Unsupported input format"**, check file suffixes.
- Tip: rename files to `.fastq(.gz)` or `.fasta(.gz)` and rerun the converter.

### The Starting Line
- If you see **"One or more rounds have zero total reads"**, normalisation cannot proceed safely.
- This means at least one round is effectively empty.
- Tip: check that all round columns are mapped correctly and not accidentally blank.

### The Race Begins
- If scores look flat (many near 0.5), trajectories may be too similar or too sparse.
- That is not a code crash, but it is a signal to inspect round quality and `min_total_count`.
- Tip: inspect `log2_enrichment`, `trend_slope`, and `terminal_guardrail` in exported results. `pace_consistency` is now kept for diagnostics rather than the main ranked output.
- If this stage is slow with very large candidate sets, try `scoring.vectorized_metrics: true`.

### Security Check
- If you see **"Failed to fetch PDB target"** or **"No sequence found"**, the run is correctly blocking unsafe analysis.
- Tip: check network access, ID spelling, or switch to a local FASTA target file.
- The run stops early here because silent fallback would corrupt later decisions.

### The Winning Bunch
- If you get **"No candidates passed filters"**, the filters are likely too strict for the current dataset.
- Tip: loosen `min_log2_enrichment` or lower `library.min_total_count` gradually and rerun.
- Keep a record of threshold changes so shortlist criteria can be justified later.
- If you see **"scoring.diversity_kmer_size must be >= 1"**, set `scoring.diversity_kmer_size` to `3` and rerun.
- If ranking still feels slow, raise `library.min_total_count` to reduce the candidate pool before Station 5.
- If ViennaRNA is not installed, nothing breaks. The optional structure annotations will simply not appear.

## Project Structure

```text
src/
├── pipeline.py            # Orchestrates the full run and stage-by-stage execution
├── sequence_generator.py  # The Scanner + The Starting Line logic
├── binding_scorer.py      # The Race Begins scoring
├── target_analyzer.py     # Security Check target validation
├── filter_rank.py         # The Winning Bunch filtering and ranking
├── structure_predictor.py # Optional structure annotations
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

It reports five universal markers:
- Leaderboard stability between the last two rounds
- Top-candidate slope trajectory (acceleration/deceleration)
- Pool dominance coverage (top 1 / 10 / 100)
- Library health in the final round (reads, unique sequences, redundancy ratio)
- Pace consistency across top-ranked candidates

Recommendation output:
- `A`: Sequence more rounds
- `B`: Stop and validate
- `C`: Potential over-selection, review earlier rounds

## License

This project is released under the MIT License. See [LICENSE](LICENSE).
