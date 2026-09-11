"""Beginner-friendly launcher for the aptamer pipeline.

This module gives us a gentle front door for teammates who do not want to type
paths or edit YAML by hand. We still save the generated config so each run stays
traceable and reproducible.
"""

from __future__ import annotations

import argparse
import logging
from copy import deepcopy
from pathlib import Path
from typing import Optional

import yaml

from src.fasta_round_counter import convert_round_files, _infer_round_number
from src.pipeline import run_pipeline
from src.utils import ensure_output_dir, load_config, setup_logging

logger = logging.getLogger("aptamer_pipeline")


def _default_base_config() -> dict:
    config_path = Path(__file__).resolve().parents[1] / "config" / "pipeline_config.yaml"
    return load_config(str(config_path))


def _prompt_choice(title: str, options: list[str]) -> int:
    print(f"\n{title}")
    for idx, option in enumerate(options, start=1):
        print(f"  {idx}. {option}")
    while True:
        raw = input("Choose a number: ").strip()
        if raw.isdigit():
            choice = int(raw)
            if 1 <= choice <= len(options):
                return choice - 1
        print("Please enter one of the listed numbers.")


def _prompt_text(prompt: str, default: Optional[str] = None) -> str:
    suffix = f" [{default}]" if default else ""
    while True:
        value = input(f"{prompt}{suffix}: ").strip()
        if value:
            return value
        if default is not None:
            return default
        print("This field is required.")


def _prompt_yes_no(prompt: str, default: bool = False) -> bool:
    default_text = "Y/n" if default else "y/N"
    while True:
        value = input(f"{prompt} [{default_text}]: ").strip().lower()
        if not value:
            return default
        if value in {"y", "yes"}:
            return True
        if value in {"n", "no"}:
            return False
        print("Please answer yes or no.")


def _file_dialog_available() -> bool:
    try:
        import tkinter  # noqa: F401
        from tkinter import filedialog  # noqa: F401
        return True
    except Exception:
        return False


def _pick_single_file(title: str) -> str:
    if _file_dialog_available():
        try:
            import tkinter as tk
            from tkinter import filedialog

            root = tk.Tk()
            root.withdraw()
            root.update()
            path = filedialog.askopenfilename(title=title)
            root.destroy()
            if path:
                return path
        except Exception:
            pass
    return _prompt_text(f"{title} (paste the full path)")


def _pick_multiple_files(title: str) -> list[str]:
    if _file_dialog_available():
        try:
            import tkinter as tk
            from tkinter import filedialog

            root = tk.Tk()
            root.withdraw()
            root.update()
            paths = filedialog.askopenfilenames(title=title)
            root.destroy()
            if paths:
                return list(paths)
        except Exception:
            pass

    print(f"{title}")
    print("Paste the full paths separated by commas.")
    raw = _prompt_text("Round files")
    return [part.strip() for part in raw.split(",") if part.strip()]


def _pick_directory(title: str, default_dir: str) -> str:
    if _file_dialog_available():
        try:
            import tkinter as tk
            from tkinter import filedialog

            root = tk.Tk()
            root.withdraw()
            root.update()
            path = filedialog.askdirectory(title=title, initialdir=default_dir)
            root.destroy()
            if path:
                return path
        except Exception:
            pass
    return _prompt_text(f"{title} (folder path)", default=default_dir)


def _build_run_config(
    base_config: dict,
    counts_file: str,
    target_input_type: str,
    target_input_value: str,
    output_dir: str,
) -> dict:
    config = deepcopy(base_config)
    config["selex"]["counts_file"] = counts_file
    config["selex"].pop("round_columns", None)
    config["target"]["input_type"] = target_input_type
    config["target"]["input_value"] = target_input_value
    config["target"]["name"] = Path(target_input_value).stem if target_input_type == "fasta" else target_input_value
    config["output"]["directory"] = output_dir
    config["output"]["generate_plots"] = False
    return config


def _write_generated_config(config: dict, output_dir: Path) -> Path:
    config_path = output_dir / "interactive_run_config.yaml"
    inputs = [Path(config["selex"]["counts_file"]).resolve()]
    if config["target"]["input_type"] == "fasta":
        inputs.append(Path(config["target"]["input_value"]).resolve())
    if config_path.resolve() in inputs:
        raise ValueError("The settings file would overwrite an input. Choose a different output folder.")
    with open(config_path, "w") as handle:
        yaml.safe_dump(config, handle, sort_keys=False)
    return config_path


def _prompt_target() -> tuple[str, str]:
    choice = _prompt_choice(
        "How would you like to provide the target?",
        [
            "Pick a target FASTA file",
            "Type a UniProt ID",
            "Type a PDB ID",
        ],
    )
    if choice == 0:
        return "fasta", _pick_single_file("Choose the target FASTA file")
    if choice == 1:
        return "uniprot", _prompt_text("Enter the UniProt ID")
    return "pdb_id", _prompt_text("Enter the PDB ID")


def _run_from_counts_mode(base_config: dict) -> None:
    counts_file = _pick_single_file("Choose the SELEX counts file (.csv or .tsv)")
    target_type, target_value = _prompt_target()
    output_dir = _pick_directory("Choose an output folder", "output/interactive_run")
    output_path = ensure_output_dir(output_dir)

    config = _build_run_config(
        base_config=base_config,
        counts_file=counts_file,
        target_input_type=target_type,
        target_input_value=target_value,
        output_dir=str(output_path),
    )
    if _prompt_yes_no("Specify a round order or use only some rounds?", default=False):
        names = _prompt_text("Round labels in selection order, separated by commas")
        config["selex"]["round_columns"] = [name.strip() for name in names.split(",")]
    config_path = _write_generated_config(config, output_path)

    print("\nStarting the pipeline from your counts table.")
    print(f"Saved run settings to: {config_path}")
    run_pipeline(config, stage="all")
    print(f"Finished. Results are in: {output_path}")


def _run_from_round_files_mode(base_config: dict) -> None:
    round_files = [Path(p) for p in _pick_multiple_files("Choose one FASTA/FASTQ file per round")]
    if len(round_files) < 2:
        raise ValueError("Please choose at least two round files.")

    target_type, target_value = _prompt_target()
    output_dir = _pick_directory("Choose an output folder", "output/interactive_run")
    output_path = ensure_output_dir(output_dir)

    left_anchor = None
    right_anchor = None
    if _prompt_yes_no("Do these reads still contain constant primer regions?", default=False):
        left_anchor = _prompt_text("Enter the left anchor sequence")
        right_anchor = _prompt_text("Enter the right anchor sequence")

    print("\nConfirm the actual selection round for each file; paired-end mates are not separate rounds.")
    round_labels = []
    for index, path in enumerate(round_files, start=1):
        inferred = _infer_round_number(path)
        default = str(inferred if inferred is not None else index)
        while True:
            number = _prompt_text(f"Round number for {path.name}", default=default)
            if number.isdigit():
                round_labels.append(f"round_{int(number)}")
                break
            print("Enter a non-negative round number, such as 0 or 3.")

    counts_file = output_path / "interactive_round_counts.csv"
    summary_file = output_path / "interactive_round_summary.tsv"

    print("\nConverting the round files into a counts table first.")
    convert_round_files(
        round_files=round_files,
        output_csv=counts_file,
        summary_tsv=summary_file,
        round_labels=round_labels,
        left_anchor=left_anchor,
        right_anchor=right_anchor,
    )

    config = _build_run_config(
        base_config=base_config,
        counts_file=str(counts_file),
        target_input_type=target_type,
        target_input_value=target_value,
        output_dir=str(output_path),
    )
    config["conversion"] = {
        "round_files": [str(path.resolve()) for path in round_files],
        "round_labels": round_labels,
        "left_anchor": left_anchor,
        "right_anchor": right_anchor,
    }
    config_path = _write_generated_config(config, output_path)

    print(f"Saved round counts to: {counts_file}")
    print(f"Saved run settings to: {config_path}")
    run_pipeline(config, stage="all")
    print(f"Finished. Results are in: {output_path}")


def main() -> None:
    parser = argparse.ArgumentParser(description="Choose files interactively, then run the usual aptamer pipeline.")
    parser.parse_args()
    setup_logging(verbose=True)
    base_config = _default_base_config()

    print("Aptamer Pipeline Interactive Launcher")
    print("We will collect the files, save a run config, and then run the usual pipeline.")

    choice = _prompt_choice(
        "What do you have right now?",
        [
            "A SELEX counts file (.csv or .tsv)",
            "One FASTA/FASTQ file per round",
        ],
    )

    if choice == 0:
        _run_from_counts_mode(base_config)
    else:
        _run_from_round_files_mode(base_config)


if __name__ == "__main__":
    main()
