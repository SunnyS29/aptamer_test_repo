import csv
import gzip

from src.fasta_round_counter import convert_fasta_rounds, convert_round_files


def test_convert_fasta_rounds_infers_rounds_and_counts(tmp_path):
    round0 = tmp_path / "selection_r0.fasta"
    round1 = tmp_path / "selection_r1.fasta"
    round2 = tmp_path / "selection_r2.fasta"

    round0.write_text(
        ">a\nACGT\n>b\nACGT\n>c\nTGCA\n"
    )
    round1.write_text(
        ">a\nACGT\n>b\nGGGG\n"
    )
    # Includes an RNA sequence (U) to confirm normalization to T.
    round2.write_text(
        ">a\nACGU\n>b\nGGGG\n>c\nGGGG\n"
    )

    out_csv = tmp_path / "counts.csv"
    out_summary = tmp_path / "summary.tsv"
    summaries = convert_fasta_rounds(
        [round2, round0, round1],  # intentionally shuffled to test inference + ordering
        output_csv=out_csv,
        round_labels=None,
        summary_tsv=out_summary,
    )

    assert [row.round_label for row in summaries] == ["round_0", "round_1", "round_2"]

    with out_csv.open() as handle:
        rows = list(csv.DictReader(handle))
    by_seq = {row["sequence"]: row for row in rows}

    assert by_seq["ACGT"]["round_0"] == "2"
    assert by_seq["ACGT"]["round_1"] == "1"
    assert by_seq["ACGT"]["round_2"] == "1"
    assert by_seq["GGGG"]["round_0"] == "0"
    assert by_seq["GGGG"]["round_1"] == "1"
    assert by_seq["GGGG"]["round_2"] == "2"

    with out_summary.open() as handle:
        summary_rows = list(csv.DictReader(handle, delimiter="\t"))
    assert len(summary_rows) == 3
    assert summary_rows[0]["round"] == "round_0"


def test_convert_fasta_rounds_with_explicit_labels(tmp_path):
    r_a = tmp_path / "alpha.fasta"
    r_b = tmp_path / "beta.fasta"
    r_a.write_text(">x\nAAAA\n>y\nCCCC\n")
    r_b.write_text(">x\nAAAA\n>x2\nAAAA\n")

    out_csv = tmp_path / "custom_counts.csv"
    convert_fasta_rounds(
        [r_a, r_b],
        output_csv=out_csv,
        round_labels=["round_start", "round_end"],
        summary_tsv=None,
    )

    with out_csv.open() as handle:
        rows = list(csv.DictReader(handle))
    assert "round_start" in rows[0]
    assert "round_end" in rows[0]
    by_seq = {row["sequence"]: row for row in rows}
    assert by_seq["AAAA"]["round_start"] == "1"
    assert by_seq["AAAA"]["round_end"] == "2"


def test_convert_round_files_from_fastq_and_fasta(tmp_path):
    round1 = tmp_path / "sample_round_1.fastq"
    round2 = tmp_path / "sample_round_2.fasta.gz"

    # FASTQ round has 3 reads, including one RNA sequence that should normalize U -> T.
    round1.write_text(
        "@r1\nACGU\n+\n####\n"
        "@r2\nGGGG\n+\n####\n"
        "@r3\nGGGG\n+\n####\n"
    )

    with gzip.open(round2, "wt") as handle:
        handle.write(">x\nACGT\n>y\nTTTT\n")

    out_csv = tmp_path / "mixed_counts.csv"
    summaries = convert_round_files(
        [round2, round1],  # intentionally shuffled; round labels should be inferred + sorted
        output_csv=out_csv,
        round_labels=None,
        summary_tsv=None,
    )

    assert [row.round_label for row in summaries] == ["round_1", "round_2"]

    with out_csv.open() as handle:
        rows = list(csv.DictReader(handle))
    by_seq = {row["sequence"]: row for row in rows}

    assert by_seq["ACGT"]["round_1"] == "1"
    assert by_seq["ACGT"]["round_2"] == "1"
    assert by_seq["GGGG"]["round_1"] == "2"
    assert by_seq["GGGG"]["round_2"] == "0"
    assert by_seq["TTTT"]["round_1"] == "0"
    assert by_seq["TTTT"]["round_2"] == "1"


def test_convert_round_files_extracts_insert_between_anchors(tmp_path):
    round1 = tmp_path / "round_1.fastq"
    round2 = tmp_path / "round_2.fastq"
    round1.write_text(
        "@r1\nLEFTAACCRIGHT\n+\n##############\n"
        "@r2\nLEFTAACCRIGHT\n+\n##############\n"
        "@r3\nLEFTGGTTRIGHT\n+\n##############\n"
        "@r4\nNOANCHOR\n+\n########\n"
    )
    round2.write_text(
        "@r1\nLEFTAACCRIGHT\n+\n##############\n"
        "@r2\nLEFTTTAARIGHT\n+\n##############\n"
    )

    out_csv = tmp_path / "anchored_counts.csv"
    summaries = convert_round_files(
        [round1, round2],
        output_csv=out_csv,
        round_labels=["round_1", "round_2"],
        summary_tsv=None,
        left_anchor="LEFT",
        right_anchor="RIGHT",
        allow_reverse_complement=False,
    )

    assert summaries[0].skipped_unmatched == 1

    with out_csv.open() as handle:
        rows = list(csv.DictReader(handle))
    by_seq = {row["sequence"]: row for row in rows}
    assert by_seq["AACC"]["round_1"] == "2"
    assert by_seq["AACC"]["round_2"] == "1"
    assert by_seq["GGTT"]["round_1"] == "1"
    assert by_seq["TTAA"]["round_2"] == "1"
