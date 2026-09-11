# Aptamer Pipeline
## An Empirical HT-SELEX Enrichment & Analysis Pipeline

This pipeline moves from raw HT-SELEX count tables to a ranked shortlist of aptamer candidates.
The pipeline ranks sequences by changes in normalized abundance across rounds rather than by abundance in a single round.
This README explains each processing step and the checks behind the final output.

## What This Pipeline Is (and Is Not)

- **It is** an empirical analysis pipeline that reads real sequence counts from SELEX rounds.
- **It is not** a random-sequence generator.
- **It records the run:** logs give filter totals, and exported scores explain the retained shortlist. We do not yet export a separate rejection reason for every discarded sequence.
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
- Uses three signals together: first-to-last log2 enrichment, overall trend slope, and a strong but graded terminal guardrail.
- Log2 enrichment and slope do the main ranking work.
- The terminal guardrail catches sequences that fade in the last round, without penalising step-function winners that take off late.

### 4. Security Check (Target Verification)
- Fetches target information from PDB/UniProt (or reads FASTA).
- If target retrieval fails, the run hard-stops.
- This checks target-file format and retrieval, not binding. Target features do not contribute to the enrichment score.

### 5. The Winning Bunch (Filtering + Ranking)
- Removes weak candidates using enrichment thresholds.
- Ranks candidates by their enrichment score after applying the configured enrichment threshold.
- Computes pooled k-mer rarity across the final shortlist as an annotation; diversity cannot change rank.
- The output shortlist is saved to CSV/JSON for downstream review.

## How The Race Works

The following sections explain what each calculation does and why we use it.
Think of every sequence as a runner in a stadium. The question is not who looked good once. It is who keeps moving toward the front as the race gets harder.

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
- Least-squares slope: fits a straight trend line across sampled rounds to measure overall upward movement. Samples are equally spaced in this calculation, even if biological rounds were skipped.
- Terminal guardrail: checks whether the last round goes down versus the round before it.
- Min-max scaling: rescales the growth features to `0-1` so they can be combined fairly.
- Weighted core: fold change `0.85` and slope `0.15`; the result is then multiplied once by the squared terminal guardrail for a strong, graded fade penalty.
- A larger final-round decline causes a larger reduction in the candidate's score.

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
- Minimum log2 enrichment: removes candidates below the configured evidence floor.
- Enrichment-score ordering: stronger measured trajectories rank first.
- Diversity score: reports k-mer rarity within the shortlist as an annotation only.
- Simple meaning: enrichment picks the winners; diversity helps us inspect whether the shortlist contains different sequence families.

### 7. The Finish-Line Referee
- Job: decide whether the experiment looks mature enough to stop.
- Exact methods used:
- Top-10 overlap percentage between the last two rounds.
- Jaccard index for set similarity of the two leaderboards.
- Coverage percentages for the top 1, top 10, and top 100 sequences.
- Mean acceleration of the top 3 trajectories.
- Mean, median, and coefficient of variation of a pace score based on monotonicity and RMSE.
- A composite data quality score and a rule-based recommendation.

## Quick Start

### Install
Use Python 3.11 for the same version checked by the automated test workflow.

```bash
git clone https://github.com/SunnyS29/aptamer_test_repo.git
cd aptamer_test_repo
pip install -r requirements.txt
```

Optional for plots (`output.generate_plots: true`):
```bash
pip install -r requirements-plots.txt
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
- let you choose an output folder and confirm the selection-round numbers for raw files
- save an `interactive_run_config.yaml` file with the chosen settings and, for raw files, the round mapping and extraction anchors

The sample config and dataset presets use a synthetic demonstration target. Replace `target.input_value` and `target.name` with your experimental target before using its exported target context. The demo header is retained in JSON output and produces a warning.

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
- If none of the file names contain round numbers, files are sorted alphabetically. Check that this is the real selection order, or provide `--round-labels`. Partly labelled or duplicate rounds now stop rather than being silently renamed.
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
- Each round has its own count column (`round_1`, `round_2`, etc.). Headers must be unique.
- Use an explicit `0` when a sequence was not observed. Blank counts and incomplete rows now stop the run, because a missing measurement is not evidence of absence.

### Option 2: Long table (one row per sequence per round)
Use this if your export is already in long format.

```csv
sequence,round,count
ACGT...,round_1,10
ACGT...,round_2,25
```

What matters:
- Must include all three columns: `sequence`, `round`, `count`.
- Counts must be whole, non-negative numbers. Present rows need a count; a missing sequence-round pair is treated as zero in this sparse format. Use it only when an omitted pair really means no reads were observed.
- Sequence text is uppercased and `U` is converted to `T`. Candidates with ambiguous or invalid bases, including `N`, are excluded by sequence QC rather than offered for synthesis.

### Option 3: Raw sequencing files (FASTQ/FASTA)
Use this when starting from raw files from the sequencer.

- Supported: `.fastq`, `.fastq.gz`, `.fasta`, `.fasta.gz`
- One file should represent one selection round, not one paired-end mate or one sequencing lane. Merge paired reads or combine lanes appropriately before using the converter.
- FASTA records must represent individual reads. A published list of winners or a file with one record per unique sequence does not preserve read counts.
- We check FASTQ record structure and matching sequence/quality lengths. We do not filter on Phred scores or correct sequencing errors; prepare reads for your experiment before counting.
- Convert first with `python -m src.fasta_round_counter ...`
- Then run the pipeline on the new counts table.

Common input mistakes:
- Missing `sequence` column name.
- A round column with all zeros, or a row with counts but no sequence. Missing sequences now stop the run rather than silently losing their counts.
- Mixed files from different experiments in one run.

## Configuration Cheat Sheet

Edit `config/pipeline_config.yaml`:

- `target`: where target info comes from (`pdb_id`, `fasta`, `smiles`, `uniprot`)
- `selex.counts_file`: the counts table path
- `selex.round_columns`: optional round names in the exact selection order to use, for either wide or long tables. For example, `[round_start, round_end]` preserves that order. Without this setting, numbered round labels are sorted numerically; custom names need an explicit order. Numeric metadata such as `length` are not treated as rounds.
- `library`: sequence QC filters (length, GC, homopolymer, min total count)
- `scoring`: pseudocount + growth weights. The pseudocount must be finite and positive; weights must be finite, non-negative, and have a positive sum.
- `scoring.vectorized_metrics`: set `true` to speed up enrichment and slope calculations with NumPy on large libraries (default `false`)
- The squared terminal guardrail is applied once after the weighted enrichment score; it is not an additional weighted component
- `scoring.diversity_kmer_size`: k-mer size used for diversity rarity scoring (default `3`)
- `filtering`: shortlist strictness. `top_n` must be a positive integer; `min_log2_enrichment` must be a finite number or `null`.
- `output`: file format + output directory

## Friendly Troubleshooting (By Section)

### The Scanner
- If the pipeline says **"Could not find a sequence column"**, it usually means column headers are inconsistent.
- Tip: rename the sequence column to `sequence` and rerun.
- If conversion fails with **"Unsupported input format"**, check file suffixes.
- Tip: check that the file really contains FASTA or FASTQ before correcting its suffix. Renaming a spreadsheet will not convert it.
- **"Round labels must be non-empty and unique"** means two inputs may have been assigned to the same round. Check the file-to-round mapping; do not label paired-end mates as different selection rounds.
- **"FASTQ sequence and quality lengths differ"** means the record is damaged or incorrectly exported. Check the source download rather than trimming characters just to make it pass.

### The Starting Line
- If you see **"One or more rounds have zero total reads"**, normalisation cannot proceed safely.
- This means at least one round is effectively empty.
- Tip: check that all round columns are mapped correctly and not accidentally blank.

### The Race Begins
- If scores look flat (many near 0.5), trajectories may be too similar or too sparse.
- That is not a code crash, but it is a signal to inspect round quality and `min_total_count`.
- Tip: inspect `log2_enrichment`, `trend_slope`, and `terminal_guardrail` in exported results. The separate stopping diagnostic calculates pace from the original round counts when needed.
- If this stage is slow with very large candidate sets, try `scoring.vectorized_metrics: true`.

### Security Check
- If you see **"Failed to fetch PDB target"** or **"No sequence found"**, the run is correctly blocking unsafe analysis.
- Tip: check network access, ID spelling, or switch to a local FASTA target file.
- If UniProt returns a web page instead of one FASTA record, the run stops even when the HTTP request says it succeeded.
- Check target identity yourself: local FASTA uses its first record, and PDB retrieval uses polymer entity 1. Valid sequence letters alone cannot confirm that you chose the intended binding protein.

### The Winning Bunch
- If you get **"No candidates passed filters"**, inspect the evidence before relaxing the filters. The completed run now writes a header-only CSV and an empty JSON result, replacing old winners and removing an old summary plot.
- Structure fields in a normal `all` run are placeholders, not measured results: `mfe=0`, `motif_count=0`, and `has_g_quadruplex=False` do not establish the absence of structure.
- Tip: loosen `min_log2_enrichment` or lower `library.min_total_count` gradually and rerun.
- Keep a record of threshold changes so shortlist criteria can be justified later.
- If you see **"scoring.diversity_kmer_size must be >= 1"**, set `scoring.diversity_kmer_size` to `3` and rerun.
- If ranking still feels slow, raise `library.min_total_count` to reduce the candidate pool before Station 5.
- If ViennaRNA is not installed, the pipeline continues without optional structure annotations.

## Run Records

A full `all` run writes `run_manifest.json` with its effective config, ordered rounds, input paths and SHA-256 fingerprints, Python version, source-code fingerprints, and Git revision when available. Hashing reads the counts file once more in small chunks, without keeping another copy in memory.

The manifest starts with `status: started` and becomes `complete` only after export and any requested plots finish. If it remains `started`, inspect the error output; that run did not finish. These fingerprints cover the counts table and a local target FASTA, not the original raw-read files.

Reusing an output folder replaces the standard ranked CSV/JSON, summary plot, default validation report, and manifest. Use a different output folder to keep an earlier run. Other files are left alone. The plots show enrichment, final CPM, and ranking score, not binding measurements or structure placeholders.

## Project Structure

```text
src/
├── interactive_launcher.py # Prompts for files, then calls the normal pipeline
├── run_record.py          # Records inputs and manages current-run output files
├── pipeline.py            # Orchestrates the full run and stage-by-stage execution
├── sequence_generator.py  # The Scanner + The Starting Line logic
├── binding_scorer.py      # The Race Begins scoring
├── target_analyzer.py     # Security Check target validation
├── filter_rank.py         # The Winning Bunch filtering and ranking
├── structure_predictor.py # Optional structure annotations
├── stopping_diagnostic.py # Optional stopping-point health check
├── validation_diagnostic.py # Optional bootstrap + walk-forward checks
└── utils.py               # Shared helpers
```

## Testing

```bash
python -m pip install -r requirements-dev.txt
python -m pytest tests/ -v
```

## Stopping-Point Diagnostic

Use this when you want a quick health check on whether SELEX rounds are converging. With only two rounds, acceleration and the composite data quality score are reported as unavailable (`null` in JSON); there is only one observed interval to compare:

```bash
python -m src.stopping_diagnostic --config config/pipeline_config.yaml
```

It reports five screening markers:
- Leaderboard stability between the last two rounds
- Top-candidate slope trajectory (acceleration/deceleration)
- Pool dominance coverage (top 1 / 10 / 100). Raw percentages use the whole count table; ranked-pool percentages use only exported shortlisted sequences.
- Library health in the final round (reads, unique sequences, redundancy ratio)
- Pace consistency across top-ranked candidates

Recommendation output:
- `A`: More evidence needed. If sampling is sparse or the leaderboard is tied, review sequencing depth before adding selection rounds.
- `B`: Stop and validate
- `C`: Potential over-selection, review earlier rounds

The stop check now reads the same CSV/TSV formats and configured rounds as the main pipeline. It excludes zero-count entries from the leaderboards and flags ties at the top-10 boundary.

The library-health route to `B` requires at least 80% overlap, fewer than 10 million observed unique sequences, and at least 5 reads per observed unique sequence on average. Low redundancy blocks a stop recommendation, and tied leaderboard boundaries cannot establish convergence. A high average can still hide many singletons or PCR duplicates; repeated reads are not independent laboratory replicates.

These cutoffs and the 0-100 data quality score are heuristics, not calibrated probabilities or universally validated lab thresholds. Dominance can prompt an over-selection warning, but it cannot establish PCR artifacts, binding affinity, or the target-to-aptamer ratio. Use the recommendation to plan a lab review, not as an automatic instruction to discard a round.

## Optional Confidence Failsafe

This is not a sixth station and it does not change **The Winning Bunch**. We run it after the normal pipeline when we want more evidence before committing laboratory time and reagents to a shortlist.

```bash
python -m src.validation_diagnostic \
  --config config/pipeline_config.yaml \
  --bootstrap-replicates 200 \
  --top-k 10
```

The report is saved as `validation_report.json` inside the configured output directory.

What it tells us:

- **Bootstrap confidence:** we resample each round at the same read depth and rerun the enrichment score. The result shows how often each original leader stays in the top `K`.
- **95% rank interval:** the middle 95% of resampled ranks shows how far a candidate moves under this sampling model. A candidate that fails the enrichment floor is assigned rank `candidate_count + 1` for that repeat; the interval is not a probability of binding.
- **Walk-forward overlap:** we hide one later round, rank candidates using only the earlier rounds, and then compare our prediction with the hidden leaders.
- **Spearman correlation:** this compares scores of candidates passing the training enrichment floor with their hidden-round CPM. `1` means strong agreement, `0` means little rank relationship, and a negative value means the order tends to reverse. With fewer than two candidates or constant values, the result is `null`, not zero, and is excluded from the correlation average. The report includes the number of candidates and usable splits.
- **Training eligibility:** this reports whether hidden leaders had enough earlier evidence to enter the race. Sparse early splits are labelled `not_evaluable` and left out of averages instead of being given a misleading zero.

Where this check stops:

- Bootstrap confidence measures read-sampling uncertainty within the original QC/count-retained candidate set. Excluded sequences stay in a pooled background, and unobserved sequences cannot appear. It cannot remove PCR bias, non-specific selection, or biological variation.
- Walk-forward validation tests future sequencing abundance, not physical target binding.
- We rebuild every training pool from its earlier rounds only, so later counts cannot leak into an earlier prediction.
- Start with `20-50` bootstrap replicates for a quick check. Use at least `200` for a final report; runtime grows roughly with the number of replicates.
- The random seed defaults to `42`, so we can reproduce the same result later. We now prepare fixed round totals and sampling probabilities once, rather than rescanning the full table for every repeat.
- Walk-forward evaluation does not undo earlier tuning on this dataset. Keep weights fixed before evaluating a new experiment if you want an independent assessment.

Helpful tips:

- **"Walk-forward validation needs at least three rounds"** means we need two rounds for scoring and one later round for validation.
- If NumPy is missing, run `pip install -r requirements.txt`.
- If no candidates pass the enrichment floor, check `filtering.min_log2_enrichment` before adding more computation.

## License

This project is released under the MIT License. See [LICENSE](LICENSE).
