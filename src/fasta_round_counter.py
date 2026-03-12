"""Turn one sequencing file per round into one count table.

This helper is the bridge between raw round files and the main pipeline.
Its job is simple:
1. open each FASTA/FASTQ file,
2. pull out the sequence we care about,
3. count how often each sequence appears,
4. write one table the rest of the pipeline can compare across rounds.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import re
from collections import Counter
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, TextIO

FASTA_SUFFIXES = (".fasta", ".fa", ".fna", ".fas")
FASTQ_SUFFIXES = (".fastq", ".fq")


@dataclass
class RoundSummary:
    round_label: str
    file_path: str
    total_reads: int
    unique_sequences: int
    skipped_empty: int = 0
    skipped_unmatched: int = 0


@dataclass
class _RoundCounts:
    """Keep one round's counts and skip totals together while we process it."""

    label: str
    file_path: str
    counter: Counter[str]
    skipped_empty: int
    skipped_unmatched: int


def _open_text(path: Path) -> TextIO:
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt")
    return path.open("r")


def _basename_without_suffixes(path: Path) -> str:
    name = path.name
    if name.lower().endswith(".gz"):
        name = name[:-3]
    for suffix in (*FASTA_SUFFIXES, *FASTQ_SUFFIXES):
        if name.lower().endswith(suffix):
            return name[: -len(suffix)]
    return Path(name).stem


def _detect_file_format(path: Path) -> str:
    name = path.name.lower()
    if name.endswith(".gz"):
        name = name[:-3]

    if name.endswith(FASTA_SUFFIXES):
        return "fasta"
    if name.endswith(FASTQ_SUFFIXES):
        return "fastq"

    raise ValueError(
        f"Unsupported input format for '{path}'. "
        "Use FASTA/FASTQ files with extensions: "
        ".fasta/.fa/.fna/.fas/.fastq/.fq (optionally .gz)."
    )


def _iter_fasta_sequences(path: Path) -> Iterable[str]:
    """Read sequences from a FASTA file.

    FASTA entries can span multiple lines, so we join them before yielding.
    """
    seq_parts: list[str] = []

    with _open_text(path) as handle:
        for line in handle:
            text = line.strip()
            if not text:
                continue
            if text.startswith(">"):
                if seq_parts:
                    yield "".join(seq_parts)
                    seq_parts = []
                continue
            seq_parts.append(text)

    if seq_parts:
        yield "".join(seq_parts)


def _iter_fastq_sequences(path: Path) -> Iterable[str]:
    """Read sequence lines from FASTQ records.

    FASTQ stores each read in blocks of four lines:
    header, sequence, plus-line, quality.
    We only need the sequence line for counting.
    """
    with _open_text(path) as handle:
        line_number = 0
        while True:
            header = handle.readline()
            if not header:
                break
            seq = handle.readline()
            plus = handle.readline()
            qual = handle.readline()
            line_number += 4

            if not (seq and plus and qual):
                raise ValueError(
                    f"Incomplete FASTQ record in '{path}' near line {line_number}."
                )
            if not header.startswith("@") or not plus.startswith("+"):
                raise ValueError(
                    f"Malformed FASTQ record in '{path}' near line {line_number}."
                )
            yield seq.strip()


def _iter_sequences(path: Path) -> Iterable[str]:
    """Pick the correct reader based on the file extension."""
    file_format = _detect_file_format(path)
    if file_format == "fasta":
        yield from _iter_fasta_sequences(path)
    else:
        yield from _iter_fastq_sequences(path)


def _normalize_sequence(seq: str) -> str:
    """Clean a sequence so small formatting differences do not split counts."""
    return re.sub(r"\s+", "", seq).upper().replace("U", "T")


def _reverse_complement(seq: str) -> str:
    """Return the reverse complement for anchor matching when reads are flipped."""
    return seq.translate(str.maketrans("ACGTN", "TGCAN"))[::-1]


def _extract_between_anchors(
    seq: str,
    left_anchor: str | None,
    right_anchor: str | None,
    allow_reverse_complement: bool,
) -> str | None:
    """Return the insert sequence we want to count for one read.

    If anchors are provided, we keep only the sequence between them.
    If no anchors are provided, we count the whole normalized read.
    """
    normalized = _normalize_sequence(seq)
    if not normalized:
        return None
    if left_anchor is None and right_anchor is None:
        return normalized
    if not left_anchor or not right_anchor:
        raise ValueError("Provide both anchor sequences or neither.")

    search_space = [normalized]
    if allow_reverse_complement:
        search_space.append(_reverse_complement(normalized))

    for candidate in search_space:
        left_index = candidate.find(left_anchor)
        if left_index < 0:
            continue
        start = left_index + len(left_anchor)
        right_index = candidate.find(right_anchor, start)
        if right_index < 0:
            continue
        insert = candidate[start:right_index]
        if insert:
            return insert
    return None


def _infer_round_number(path: Path) -> int | None:
    """Try to pull a round number from common file-name patterns."""
    name = _basename_without_suffixes(path)
    patterns = [
        r"(?:^|[_\-])(?:round|rnd|r)[_\-]?(\d+)(?:$|[_\-])",
        r"(?:^|[_\-])(\d+)[rR](?:$|[_\-])",
    ]
    for pattern in patterns:
        match = re.search(pattern, name, flags=re.IGNORECASE)
        if match:
            return int(match.group(1))
    return None


def _resolve_round_labels(
    round_files: list[Path], round_labels: list[str] | None
) -> tuple[list[Path], list[str]]:
    """Decide the order and names of rounds.

    We prefer explicit labels from the user. If they are missing, we try to
    infer round numbers from file names. If that fails, we fall back to
    alphabetical order so the tool still behaves predictably.
    """
    if round_labels is not None:
        if len(round_labels) != len(round_files):
            raise ValueError(
                "When provided, --round-labels must have the same length as input files."
            )
        return round_files, round_labels

    inferred = [_infer_round_number(path) for path in round_files]
    all_inferred = all(value is not None for value in inferred)
    unique_inferred = len(set(inferred)) == len(inferred)
    if all_inferred and unique_inferred:
        ordered = sorted(
            zip(round_files, inferred), key=lambda item: int(item[1])  # type: ignore[arg-type]
        )
        ordered_files = [item[0] for item in ordered]
        ordered_labels = [f"round_{int(item[1])}" for item in ordered]
        return ordered_files, ordered_labels

    ordered_files = sorted(round_files)
    ordered_labels = [f"round_{idx + 1}" for idx in range(len(ordered_files))]
    return ordered_files, ordered_labels


def _validate_round_inputs(round_files: list[Path]) -> None:
    """Fail early on missing files or too few rounds."""
    if len(round_files) < 2:
        raise ValueError("Provide at least two round files so round comparison is meaningful.")

    missing = [str(path) for path in round_files if not path.exists()]
    if missing:
        raise FileNotFoundError("Missing input file(s): " + ", ".join(missing))


def _normalized_anchors(
    left_anchor: str | None, right_anchor: str | None
) -> tuple[str | None, str | None]:
    """Normalize anchor text once so we do not repeat that work per read."""
    normalized_left = _normalize_sequence(left_anchor) if left_anchor is not None else None
    normalized_right = _normalize_sequence(right_anchor) if right_anchor is not None else None
    return normalized_left, normalized_right


def _count_round_sequences(
    path: Path,
    label: str,
    left_anchor: str | None,
    right_anchor: str | None,
    allow_reverse_complement: bool,
) -> _RoundCounts:
    """Read one round file and count the inserts we want to keep."""
    counter: Counter[str] = Counter()
    skipped_empty = 0
    skipped_unmatched = 0

    for raw_seq in _iter_sequences(path):
        seq = _extract_between_anchors(
            raw_seq,
            left_anchor=left_anchor,
            right_anchor=right_anchor,
            allow_reverse_complement=allow_reverse_complement,
        )
        if seq is None:
            normalized_raw = _normalize_sequence(raw_seq)
            if normalized_raw:
                skipped_unmatched += 1
            else:
                skipped_empty += 1
            continue
        if not seq:
            skipped_empty += 1
            continue
        counter[seq] += 1

    return _RoundCounts(
        label=label,
        file_path=str(path),
        counter=counter,
        skipped_empty=skipped_empty,
        skipped_unmatched=skipped_unmatched,
    )


def _build_round_summaries(round_counts: list[_RoundCounts]) -> list[RoundSummary]:
    """Convert internal per-round counters into the simpler export summary."""
    return [
        RoundSummary(
            round_label=round_count.label,
            file_path=round_count.file_path,
            total_reads=sum(round_count.counter.values()),
            unique_sequences=len(round_count.counter),
            skipped_empty=round_count.skipped_empty,
            skipped_unmatched=round_count.skipped_unmatched,
        )
        for round_count in round_counts
    ]


def _sorted_sequences_for_export(
    counts_by_round: dict[str, Counter[str]], labels: list[str]
) -> list[str]:
    """Sort sequences so the most useful rows appear first in the output table.

    We rank by final-round count first, then total count across all rounds.
    That makes the exported file easier to inspect by eye.
    """
    all_sequences: set[str] = set()
    for counter in counts_by_round.values():
        all_sequences.update(counter.keys())

    totals_by_seq = {
        seq: sum(counts_by_round[label].get(seq, 0) for label in labels) for seq in all_sequences
    }
    final_label = labels[-1]
    return sorted(
        all_sequences,
        key=lambda seq: (
            -counts_by_round[final_label].get(seq, 0),
            -totals_by_seq[seq],
            seq,
        ),
    )


def _write_counts_table(
    output_csv: Path,
    counts_by_round: dict[str, Counter[str]],
    labels: list[str],
    sorted_sequences: list[str],
) -> None:
    """Write the wide count table used by the rest of the pipeline."""
    output_csv.parent.mkdir(parents=True, exist_ok=True)
    with output_csv.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["sequence", *labels])
        for seq in sorted_sequences:
            writer.writerow([seq, *[counts_by_round[label].get(seq, 0) for label in labels]])


def _write_summary_table(summary_tsv: Path, summaries: list[RoundSummary]) -> None:
    """Write a small per-round report for quick troubleshooting."""
    summary_tsv.parent.mkdir(parents=True, exist_ok=True)
    with summary_tsv.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(
            [
                "round",
                "file",
                "total_reads",
                "unique_sequences",
                "skipped_empty_sequences",
                "skipped_unmatched_sequences",
            ]
        )
        for row in summaries:
            writer.writerow(
                [
                    row.round_label,
                    row.file_path,
                    row.total_reads,
                    row.unique_sequences,
                    row.skipped_empty,
                    row.skipped_unmatched,
                ]
            )


def convert_round_files(
    round_files: list[Path],
    output_csv: Path,
    round_labels: list[str] | None = None,
    summary_tsv: Path | None = None,
    left_anchor: str | None = None,
    right_anchor: str | None = None,
    allow_reverse_complement: bool = True,
) -> list[RoundSummary]:
    """Convert round files into one sequence-by-round count table.

    This is the main entry point used by both the CLI and the interactive
    launcher. The helpers above keep each part small so it is easier to follow:
    validate inputs, count each round, write the matrix, then write a summary.
    """
    _validate_round_inputs(round_files)
    ordered_files, labels = _resolve_round_labels(round_files, round_labels)
    normalized_left, normalized_right = _normalized_anchors(left_anchor, right_anchor)

    round_counts = [
        _count_round_sequences(
            path=path,
            label=label,
            left_anchor=normalized_left,
            right_anchor=normalized_right,
            allow_reverse_complement=allow_reverse_complement,
        )
        for path, label in zip(ordered_files, labels)
    ]
    counts_by_round = {round_count.label: round_count.counter for round_count in round_counts}
    summaries = _build_round_summaries(round_counts)
    sorted_sequences = _sorted_sequences_for_export(counts_by_round, labels)
    _write_counts_table(output_csv, counts_by_round, labels, sorted_sequences)

    if summary_tsv is not None:
        _write_summary_table(summary_tsv, summaries)

    return summaries


def convert_fasta_rounds(
    fasta_files: list[Path],
    output_csv: Path,
    round_labels: list[str] | None = None,
    summary_tsv: Path | None = None,
    left_anchor: str | None = None,
    right_anchor: str | None = None,
    allow_reverse_complement: bool = True,
) -> list[RoundSummary]:
    """Backward-compatible wrapper kept for older imports and tests."""
    return convert_round_files(
        round_files=fasta_files,
        output_csv=output_csv,
        round_labels=round_labels,
        summary_tsv=summary_tsv,
        left_anchor=left_anchor,
        right_anchor=right_anchor,
        allow_reverse_complement=allow_reverse_complement,
    )


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Convert per-round FASTA/FASTQ files into sequence x round count table."
    )
    parser.add_argument(
        "round_files",
        nargs="+",
        help=(
            "Input round files (one file per round). "
            "Supported: FASTA/FASTQ with optional .gz compression."
        ),
    )
    parser.add_argument(
        "--output",
        "-o",
        required=True,
        help="Output CSV path (wide format: sequence, round_1, round_2, ...).",
    )
    parser.add_argument(
        "--round-labels",
        nargs="+",
        default=None,
        help="Optional custom round labels matching file order (e.g., round_0 round_1 round_2).",
    )
    parser.add_argument(
        "--summary",
        default=None,
        help="Optional TSV path with per-round totals and file mapping.",
    )
    parser.add_argument(
        "--left-anchor",
        default=None,
        help="Optional constant sequence immediately before the variable region.",
    )
    parser.add_argument(
        "--right-anchor",
        default=None,
        help="Optional constant sequence immediately after the variable region.",
    )
    parser.add_argument(
        "--no-revcomp",
        action="store_true",
        help="Disable reverse-complement anchor matching.",
    )
    return parser


def main() -> None:
    parser = build_arg_parser()
    args = parser.parse_args()

    round_files = [Path(path) for path in args.round_files]
    output_csv = Path(args.output)
    summary_tsv = Path(args.summary) if args.summary else None

    summaries = convert_round_files(
        round_files=round_files,
        output_csv=output_csv,
        round_labels=args.round_labels,
        summary_tsv=summary_tsv,
        left_anchor=args.left_anchor,
        right_anchor=args.right_anchor,
        allow_reverse_complement=not args.no_revcomp,
    )

    print(f"Wrote counts table: {output_csv}")
    print("Rounds processed:")
    for row in summaries:
        print(
            f"  {row.round_label}: total_reads={row.total_reads}, "
            f"unique_sequences={row.unique_sequences}, "
            f"skipped_unmatched={row.skipped_unmatched}, file={row.file_path}"
        )
    if summary_tsv:
        print(f"Wrote summary table: {summary_tsv}")


if __name__ == "__main__":
    main()
