"""Convert per-round FASTA/FASTQ files into a sequence x round count table.

We use this utility when collaborators provide one file per SELEX round and we
need a counts matrix compatible with The Scanner/The Starting Line.
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
    """Yield sequences from FASTA (supports multiline entries)."""
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
    """Yield sequence lines from FASTQ records (4-line chunks)."""
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
    file_format = _detect_file_format(path)
    if file_format == "fasta":
        yield from _iter_fasta_sequences(path)
    else:
        yield from _iter_fastq_sequences(path)


def _normalize_sequence(seq: str) -> str:
    return re.sub(r"\s+", "", seq).upper().replace("U", "T")


def _reverse_complement(seq: str) -> str:
    return seq.translate(str.maketrans("ACGTN", "TGCAN"))[::-1]


def _extract_between_anchors(
    seq: str,
    left_anchor: str | None,
    right_anchor: str | None,
    allow_reverse_complement: bool,
) -> str | None:
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


def convert_round_files(
    round_files: list[Path],
    output_csv: Path,
    round_labels: list[str] | None = None,
    summary_tsv: Path | None = None,
    left_anchor: str | None = None,
    right_anchor: str | None = None,
    allow_reverse_complement: bool = True,
) -> list[RoundSummary]:
    if len(round_files) < 2:
        raise ValueError("Provide at least two round files so round comparison is meaningful.")

    missing = [str(path) for path in round_files if not path.exists()]
    if missing:
        raise FileNotFoundError("Missing input file(s): " + ", ".join(missing))

    ordered_files, labels = _resolve_round_labels(round_files, round_labels)
    counts_by_round: dict[str, Counter[str]] = {}
    summaries: list[RoundSummary] = []
    normalized_left = _normalize_sequence(left_anchor) if left_anchor is not None else None
    normalized_right = _normalize_sequence(right_anchor) if right_anchor is not None else None

    for path, label in zip(ordered_files, labels):
        counter: Counter[str] = Counter()
        skipped_empty = 0
        skipped_unmatched = 0

        for raw_seq in _iter_sequences(path):
            seq = _extract_between_anchors(
                raw_seq,
                left_anchor=normalized_left,
                right_anchor=normalized_right,
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

        counts_by_round[label] = counter
        summaries.append(
            RoundSummary(
                round_label=label,
                file_path=str(path),
                total_reads=sum(counter.values()),
                unique_sequences=len(counter),
                skipped_empty=skipped_empty,
                skipped_unmatched=skipped_unmatched,
            )
        )

    all_sequences: set[str] = set()
    for counter in counts_by_round.values():
        all_sequences.update(counter.keys())

    output_csv.parent.mkdir(parents=True, exist_ok=True)
    totals_by_seq = {
        seq: sum(counts_by_round[label].get(seq, 0) for label in labels) for seq in all_sequences
    }
    final_label = labels[-1]
    sorted_sequences = sorted(
        all_sequences,
        key=lambda seq: (
            -counts_by_round[final_label].get(seq, 0),
            -totals_by_seq[seq],
            seq,
        ),
    )

    with output_csv.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["sequence", *labels])
        for seq in sorted_sequences:
            writer.writerow([seq, *[counts_by_round[label].get(seq, 0) for label in labels]])

    if summary_tsv is not None:
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
    """Backward-compatible wrapper for existing imports/tests."""
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
