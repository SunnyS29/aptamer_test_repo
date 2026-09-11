"""Station 1 + Station 2: The Scanner and The Starting Line.

This module does two jobs for us:
1) The Scanner: read real HT-SELEX count files and build a clean per-sequence table.
2) The Starting Line: normalize each round to CPM so comparisons are depth-aware.
"""

import csv
import logging
import re
from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

from src.utils import gc_content, has_homopolymer

logger = logging.getLogger("aptamer_pipeline")


@dataclass
class AptamerCandidate:
    """Represents an observed aptamer candidate from SELEX sequencing."""
    id: str
    sequence: str
    length: int
    gc: float
    round_counts: dict[str, int] = field(default_factory=dict)
    round_cpm: dict[str, float] = field(default_factory=dict)
    round_order: list[str] = field(default_factory=list)

    def to_dict(self) -> dict:
        return {
            "id": self.id,
            "sequence": self.sequence,
            "length": self.length,
            "gc_content": round(self.gc, 4),
            "round_counts": self.round_counts,
            "round_cpm": {k: round(v, 3) for k, v in self.round_cpm.items()},
        }


def validate_sequence(seq: str, gc_min: float, gc_max: float,
                      max_homopolymer: int) -> bool:
    """Check whether a sequence passes basic lab-friendly quality filters."""
    # Ambiguous bases cannot identify an exact candidate for synthesis.
    if not seq or any(base not in "ACGT" for base in seq):
        return False
    gc = gc_content(seq)
    if gc < gc_min or gc > gc_max:
        return False
    if has_homopolymer(seq, max_homopolymer):
        return False
    return True


def _normalize_sequence(value: str) -> str:
    """Normalize sequence text so we compare sequences consistently across files."""
    sequence = value.upper().replace("U", "T").replace(" ", "").strip()
    if not sequence:
        raise ValueError(
            "Found a counts row with no sequence. "
            "Tip: restore its sequence from the source file; dropping the row changes the round totals."
        )
    return sequence


def _parse_nonnegative_int(value: str, field_name: str) -> int:
    """Parse count fields safely.

    We require integer-like, non-negative counts because fractional or negative values
    almost always indicate a parsing issue upstream.
    """
    if value is None or str(value).strip() == "":
        raise ValueError(
            f"Field '{field_name}' has a missing count. "
            "Tip: use 0 for a confirmed absence; check the source file for missing measurements."
        )
    text = str(value).strip()

    try:
        numeric = float(text)
    except ValueError as exc:
        raise ValueError(f"Field '{field_name}' contains non-numeric value '{value}'.") from exc

    if numeric < 0:
        raise ValueError(f"Field '{field_name}' contains negative count '{value}'.")
    if not numeric.is_integer():
        raise ValueError(f"Field '{field_name}' must be integer-like, got '{value}'.")
    return int(numeric)


def _select_rounds(available: list[str], requested: Optional[list[str]] = None) -> list[str]:
    """Use the requested order, or sort round numbers when the order is clear."""
    if requested is not None:
        if not requested or any(
            not isinstance(name, str) or not name.strip() for name in requested
        ):
            raise ValueError("selex.round_columns must contain non-empty round names.")
        if len(set(requested)) != len(requested):
            raise ValueError("selex.round_columns contains duplicate rounds.")
        missing = [name for name in requested if name not in available]
        if missing:
            raise ValueError(f"Configured rounds are missing from the counts file: {', '.join(missing)}")
        return list(requested)

    numbered = []
    for name in available:
        match = re.search(r"(\d+)", name)
        if match is None:
            raise ValueError(
                "Round order is unclear. Tip: list the round names in selection order "
                "under selex.round_columns, or use round_1, round_2, and so on."
            )
        numbered.append((int(match.group(1)), name))
    if len({number for number, name in numbered}) != len(numbered):
        raise ValueError("Round numbers are duplicated. Tip: check the labels or set selex.round_columns explicitly.")
    return [name for number, name in sorted(numbered)]


def _read_counts_rows(path: str) -> tuple[list[str], list[dict[str, str]]]:
    """Read CSV/TSV count rows.

    We sniff the delimiter so collaborators can hand us either comma- or tab-separated
    exports without changing code.
    """
    filepath = Path(path)
    if not filepath.exists():
        raise FileNotFoundError(
            f"SELEX counts file not found: {path}. "
            "Tip: check 'selex.counts_file' and run from the project root."
        )

    try:
        with open(filepath, newline="") as handle:
            sample = handle.read(4096)
            handle.seek(0)

            try:
                dialect = csv.Sniffer().sniff(sample, delimiters=",\t;")
            except csv.Error:
                dialect = csv.excel

            reader = csv.DictReader(handle, dialect=dialect)
            fieldnames = [name.strip() for name in (reader.fieldnames or [])]
            if any(not name for name in fieldnames) or len({name.lower() for name in fieldnames}) != len(fieldnames):
                raise ValueError("Counts file has empty or duplicate headers. Tip: give each column a unique name.")
            reader.fieldnames = fieldnames
            rows = []
            for row in reader:
                if row is None:
                    continue
                if None in row or any(value is None for value in row.values()):
                    raise ValueError(
                        f"Counts row near line {reader.line_num} does not match the headers. "
                        "Tip: check for missing counts or an extra separator."
                    )
                cleaned = {k: (v.strip() if isinstance(v, str) else "") for k, v in row.items()}
                if any(cleaned.values()):
                    rows.append(cleaned)

            if not fieldnames:
                raise ValueError("Counts file is missing a header row.")
            if not rows:
                raise ValueError("Counts file has no data rows.")
            return fieldnames, rows
    except OSError as exc:
        raise ValueError(f"Failed to read SELEX counts file '{path}': {exc}") from exc


def _detect_sequence_column(fieldnames: list[str]) -> str:
    """Detect sequence column name from common conventions."""
    lower_to_original = {c.lower(): c for c in fieldnames}
    for key in ("sequence", "seq", "aptamer_sequence"):
        if key in lower_to_original:
            return lower_to_original[key]
    raise ValueError(
        "Could not find a sequence column. Expected one of: "
        "'sequence', 'seq', 'aptamer_sequence'. "
        "Tip: rename your sequence header to 'sequence'."
    )


def _detect_round_columns(
    fieldnames: list[str],
    round_prefix: Optional[str],
    round_columns: Optional[list[str]],
) -> list[str]:
    """Recognize round labels, never unrelated numeric metadata such as length."""
    if round_columns is not None:
        return _select_rounds(fieldnames, round_columns)
    if round_prefix:
        matching = [name for name in fieldnames if name.lower().startswith(round_prefix.lower())]
        if matching:
            return _select_rounds(matching)
    matching = [
        name for name in fieldnames
        if re.fullmatch(r"(?:round|rnd|r)[_-]?\d+", name, re.I)
    ]
    if matching:
        return _select_rounds(matching)
    raise ValueError(
        "Could not infer SELEX round columns. Tip: provide selex.round_columns "
        "in selection order, or use headers such as round_1 and round_2."
    )


def _prepare_long_format(
    rows: list[dict[str, str]],
    sequence_col: str,
    round_col: str,
    count_col: str,
    round_columns: Optional[list[str]] = None,
) -> tuple[list[dict[str, int]], list[str]]:
    """Convert long-format rows into a sequence-by-round matrix.

    Duplicate rows are summed so technical replicate merges do not drop signal.
    """
    seq_round_counts: dict[str, dict[str, int]] = defaultdict(lambda: defaultdict(int))
    seen_rounds = set()

    for row in rows:
        sequence = _normalize_sequence(row.get(sequence_col, ""))

        round_name = str(row.get(round_col, "")).strip()
        if not round_name:
            raise ValueError("Found empty round label in long-format counts table.")

        count = _parse_nonnegative_int(row.get(count_col, "0"), count_col)
        seq_round_counts[sequence][round_name] += count
        seen_rounds.add(round_name)

    if not seq_round_counts:
        raise ValueError("No valid sequences found in long-format counts table.")

    rounds = _select_rounds(list(seen_rounds), round_columns)
    table = []
    for seq, count_map in seq_round_counts.items():
        record = {"sequence": seq}
        for round_name in rounds:
            record[round_name] = int(count_map.get(round_name, 0))
        table.append(record)

    return table, rounds


def _prepare_wide_format(
    rows: list[dict[str, str]],
    sequence_col: str,
    round_prefix: Optional[str],
    round_columns: Optional[list[str]],
    fieldnames: list[str],
) -> tuple[list[dict[str, int]], list[str]]:
    """Validate wide-format rows and consolidate duplicates by sequence."""
    rounds = _detect_round_columns(fieldnames, round_prefix, round_columns)
    seq_round_counts: dict[str, dict[str, int]] = defaultdict(lambda: defaultdict(int))

    for row in rows:
        sequence = _normalize_sequence(row.get(sequence_col, ""))

        for round_name in rounds:
            count = _parse_nonnegative_int(row.get(round_name, "0"), round_name)
            seq_round_counts[sequence][round_name] += count

    if not seq_round_counts:
        raise ValueError("No valid sequences found in wide-format counts table.")

    table = []
    for seq, count_map in seq_round_counts.items():
        record = {"sequence": seq}
        for round_name in rounds:
            record[round_name] = int(count_map.get(round_name, 0))
        table.append(record)

    return table, rounds


def load_selex_counts(config: dict) -> tuple[list[dict[str, int]], list[str]]:
    """Load and validate SELEX counts before scoring.

    We require at least two rounds and non-zero round totals, otherwise enrichment
    claims would be mathematically weak or misleading.
    """
    selex_config = config.get("selex", {})
    counts_file = selex_config.get("counts_file")
    if not counts_file:
        raise ValueError(
            "SELEX counts file is required. Set 'selex.counts_file' in config."
        )

    round_prefix = selex_config.get("round_prefix", "round_")
    round_columns = selex_config.get("round_columns")
    if round_columns is not None and not isinstance(round_columns, list):
        raise ValueError("selex.round_columns must be a list when provided.")

    fieldnames, rows = _read_counts_rows(counts_file)
    sequence_col = _detect_sequence_column(fieldnames)
    lower_cols = {c.lower(): c for c in fieldnames}

    if "round" in lower_cols and "count" in lower_cols:
        table, rounds = _prepare_long_format(
            rows,
            sequence_col=sequence_col,
            round_col=lower_cols["round"],
            count_col=lower_cols["count"],
            round_columns=round_columns,
        )
    else:
        table, rounds = _prepare_wide_format(
            rows,
            sequence_col=sequence_col,
            round_prefix=round_prefix,
            round_columns=round_columns,
            fieldnames=fieldnames,
        )

    if len(rounds) < 2:
        raise ValueError(
            "At least two SELEX rounds are required for enrichment scoring. "
            "Tip: include at least two round columns (e.g., round_1 and round_2)."
        )

    totals = {r: sum(record[r] for record in table) for r in rounds}
    empty_rounds = [r for r, total in totals.items() if total <= 0]
    if empty_rounds:
        raise ValueError(
            "One or more rounds have zero total reads: "
            f"{', '.join(empty_rounds)}. Cannot normalize CPM. "
            "Tip: check round-column mapping and empty columns."
        )

    logger.info(
        "Loaded SELEX count table: %d unique sequences across %d rounds.",
        len(table), len(rounds)
    )
    logger.info(
        "Round totals: %s",
        ", ".join(f"{r}={totals[r]}" for r in rounds),
    )

    return table, rounds


def build_candidates_from_counts(
    table: list[dict], rounds: list[str], config: dict
) -> list[AptamerCandidate]:
    """Apply the pipeline's QC and CPM rules to a selected set of rounds.

    Keeping this logic in one place matters for walk-forward validation: each
    historical split must rebuild its candidate pool without looking at future
    rounds, while using exactly the same rules as a normal pipeline run.
    """
    if len(rounds) < 2:
        raise ValueError("At least two rounds are required to build candidates.")

    lib_config = config["library"]
    length_min = lib_config.get("length_min", 0)
    length_max = lib_config.get("length_max", 10_000)
    gc_min = lib_config.get("gc_min", 0.0)
    gc_max = lib_config.get("gc_max", 1.0)
    max_homo = lib_config.get("max_homopolymer", 1000)
    min_total_count = lib_config.get("min_total_count", 1)

    round_totals = {r: float(sum(record[r] for record in table)) for r in rounds}
    empty_rounds = [r for r, total in round_totals.items() if total <= 0]
    if empty_rounds:
        raise ValueError(
            "One or more selected rounds have zero total reads: "
            f"{', '.join(empty_rounds)}. Cannot normalize CPM."
        )

    logger.info(
        "Constructing candidates from observed sequences with QC filters "
        "(length=%s-%s, GC=%.2f-%.2f, min_total_count=%d).",
        length_min, length_max, gc_min, gc_max, min_total_count
    )

    candidates = []
    filtered_qc = 0
    filtered_count = 0
    end_round = rounds[-1]

    for row in table:
        seq = row["sequence"]
        length = len(seq)
        if length < length_min or length > length_max:
            filtered_qc += 1
            continue

        if not validate_sequence(seq, gc_min, gc_max, max_homo):
            filtered_qc += 1
            continue

        counts = {r: int(row[r]) for r in rounds}
        total_count = sum(counts.values())
        if total_count < min_total_count:
            filtered_count += 1
            continue

        cpm = {r: (counts[r] / round_totals[r]) * 1e6 for r in rounds}
        candidates.append(AptamerCandidate(
            id="",
            sequence=seq,
            length=length,
            gc=gc_content(seq),
            round_counts=counts,
            round_cpm=cpm,
            round_order=list(rounds),
        ))

    if not candidates:
        raise ValueError(
            "No sequences remain after SELEX ingestion and QC filters. "
            "Tip: relax QC filters or inspect your counts file for malformed rows."
        )

    candidates.sort(
        key=lambda c: (
            c.round_counts.get(end_round, 0),
            sum(c.round_counts.values()),
        ),
        reverse=True
    )
    for i, candidate in enumerate(candidates, start=1):
        candidate.id = f"APT_{i:06d}"

    logger.info(
        "Candidate construction complete: %d retained, %d dropped by QC, %d dropped by "
        "min_total_count.",
        len(candidates), filtered_qc, filtered_count
    )
    return candidates


def generate_library(config: dict) -> list[AptamerCandidate]:
    """Build observed candidates from every round in the configured count file."""
    table, rounds = load_selex_counts(config)
    return build_candidates_from_counts(table, rounds, config)
