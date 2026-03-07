"""Station 1 + Station 2: The Scanner and The Race Begins.

This module does two jobs for us:
1) The Scanner: read real HT-SELEX count files and build a clean per-sequence table.
2) The Race Begins: normalize each round to CPM so comparisons are depth-aware.
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
    gc = gc_content(seq)
    if gc < gc_min or gc > gc_max:
        return False
    if has_homopolymer(seq, max_homopolymer):
        return False
    return True


def _normalize_sequence(value: str) -> str:
    """Normalize sequence text so we compare sequences consistently across files."""
    return value.upper().replace("U", "T").replace(" ", "").strip()


def _parse_nonnegative_int(value: str, field_name: str) -> int:
    """Parse count fields safely.

    We require integer-like, non-negative counts because fractional or negative values
    almost always indicate a parsing issue upstream.
    """
    if value is None:
        return 0

    text = str(value).strip()
    if text == "":
        return 0

    try:
        numeric = float(text)
    except ValueError as exc:
        raise ValueError(f"Field '{field_name}' contains non-numeric value '{value}'.") from exc

    if numeric < 0:
        raise ValueError(f"Field '{field_name}' contains negative count '{value}'.")
    if not numeric.is_integer():
        raise ValueError(f"Field '{field_name}' must be integer-like, got '{value}'.")
    return int(numeric)


def _sort_round_names(round_names: list[str]) -> list[str]:
    """Sort round labels naturally (R1, round_2, R10) to preserve time order."""

    def sort_key(name: str):
        match = re.search(r"(\d+)", name)
        if match:
            return (0, int(match.group(1)), name)
        return (1, float("inf"), name)

    return sorted(round_names, key=sort_key)


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
            fieldnames = reader.fieldnames or []
            rows = []
            for row in reader:
                if row is None:
                    continue
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


def _is_numeric_column(rows: list[dict[str, str]], column: str) -> bool:
    """Check whether a column contains numeric-like values."""
    has_non_empty = False
    for row in rows:
        value = row.get(column, "")
        if str(value).strip() == "":
            continue
        has_non_empty = True
        try:
            float(value)
        except ValueError:
            return False
    return has_non_empty


def _detect_round_columns(
    fieldnames: list[str],
    rows: list[dict[str, str]],
    sequence_col: str,
    round_prefix: Optional[str],
    round_columns: Optional[list[str]],
) -> list[str]:
    """Detect round columns in wide-format tables.

    We try explicit config first, then prefix matching, then numeric fallback.
    This keeps ingestion flexible while still failing loudly when detection is ambiguous.
    """
    if round_columns:
        missing = [c for c in round_columns if c not in fieldnames]
        if missing:
            raise ValueError(
                "Configured round columns are missing from counts file: "
                f"{', '.join(missing)}"
            )
        return _sort_round_names(round_columns)

    if round_prefix:
        prefix_cols = [c for c in fieldnames if c.lower().startswith(round_prefix.lower())]
        if prefix_cols:
            return _sort_round_names(prefix_cols)

    excluded = {sequence_col.lower(), "aptamer_id", "id", "name", "round", "count"}
    numeric_cols = [
        col for col in fieldnames
        if col.lower() not in excluded and _is_numeric_column(rows, col)
    ]
    if len(numeric_cols) >= 2:
        return _sort_round_names(numeric_cols)

    raise ValueError(
        "Could not infer SELEX round columns. Provide 'selex.round_columns' "
        "or use a 'round_' prefix. Tip: rename columns like round_1, round_2, round_3."
    )


def _prepare_long_format(
    rows: list[dict[str, str]],
    sequence_col: str,
    round_col: str,
    count_col: str,
) -> tuple[list[dict[str, int]], list[str]]:
    """Convert long-format rows into a sequence-by-round matrix.

    Duplicate rows are summed so technical replicate merges do not drop signal.
    """
    seq_round_counts: dict[str, dict[str, int]] = defaultdict(lambda: defaultdict(int))
    seen_rounds = set()

    for row in rows:
        sequence = _normalize_sequence(row.get(sequence_col, ""))
        if not sequence:
            continue

        round_name = str(row.get(round_col, "")).strip()
        if not round_name:
            raise ValueError("Found empty round label in long-format counts table.")

        count = _parse_nonnegative_int(row.get(count_col, "0"), count_col)
        seq_round_counts[sequence][round_name] += count
        seen_rounds.add(round_name)

    if not seq_round_counts:
        raise ValueError("No valid sequences found in long-format counts table.")

    rounds = _sort_round_names(list(seen_rounds))
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
    rounds = _detect_round_columns(fieldnames, rows, sequence_col, round_prefix, round_columns)
    seq_round_counts: dict[str, dict[str, int]] = defaultdict(lambda: defaultdict(int))

    for row in rows:
        sequence = _normalize_sequence(row.get(sequence_col, ""))
        if not sequence:
            continue

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


def _load_selex_counts(config: dict) -> tuple[list[dict[str, int]], list[str]]:
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


def generate_library(config: dict) -> list[AptamerCandidate]:
    """Build observed aptamer candidates from real SELEX counts.

    Args:
        config: Full pipeline config.

    Returns:
        List of AptamerCandidate objects that pass quality filters.
    """
    lib_config = config["library"]
    length_min = lib_config.get("length_min", 0)
    length_max = lib_config.get("length_max", 10_000)
    gc_min = lib_config.get("gc_min", 0.0)
    gc_max = lib_config.get("gc_max", 1.0)
    max_homo = lib_config.get("max_homopolymer", 1000)
    min_total_count = lib_config.get("min_total_count", 1)

    table, rounds = _load_selex_counts(config)
    # We keep round totals around to convert raw counts to CPM (The Race Begins station).
    round_totals = {r: float(sum(record[r] for record in table)) for r in rounds}

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

        # CPM lets us compare rounds with different sequencing depths on equal footing.
        cpm = {r: (counts[r] / round_totals[r]) * 1e6 for r in rounds}

        candidates.append(AptamerCandidate(
            id="",
            sequence=seq,
            length=length,
            gc=gc_content(seq),
            round_counts=counts,
            round_cpm=cpm,
            round_order=rounds,
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
