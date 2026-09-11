"""Aptamer Target Identification Pipeline — CLI orchestrator.

Think of this file as the route map through our five stations:
1) The Scanner
2) The Starting Line
3) The Race Begins
4) Security Check
5) The Winning Bunch

Usage:
    python -m src.pipeline --config config/pipeline_config.yaml
    python -m src.pipeline --config config/pipeline_config.yaml --stage library
"""

import argparse
import csv
import json
import logging
import sys
from dataclasses import fields
from pathlib import Path

from src.utils import load_config, setup_logging, ensure_output_dir
from src.sequence_generator import generate_library
from src.target_analyzer import analyze_target
from src.structure_predictor import predict_structures
from src.binding_scorer import score_binding
from src.filter_rank import RankedCandidate, filter_and_rank
from src.run_record import start_run_record, finish_run_record

logger = logging.getLogger("aptamer_pipeline")

STAGES = ["target", "library", "structure", "scoring", "filtering", "all"]


def generate_plots(ranked_candidates: list, output_dir: Path) -> None:
    """Plot measured enrichment, abundance, and the score that determines rank."""
    if not ranked_candidates:
        return
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError as exc:
        logger.warning("Plots skipped: %s. Tip: install requirements-plots.txt or leave plotting off.", exc)
        return

    fig, axes = plt.subplots(1, 3, figsize=(15, 5))
    scores = [candidate.composite_score for candidate in ranked_candidates]
    enrichment = [candidate.log2_enrichment for candidate in ranked_candidates]
    abundance = [candidate.final_round_cpm for candidate in ranked_candidates]
    axes[0].hist(scores, bins=min(20, len(scores)), color="steelblue")
    axes[0].set(title="Enrichment ranking scores", xlabel="Score (not binding probability)", ylabel="Candidates")
    axes[1].scatter(enrichment, abundance, color="teal", alpha=0.7)
    axes[1].set(title="Enrichment and final abundance", xlabel="First-to-last log2 enrichment", ylabel="Final-round CPM")
    top10 = ranked_candidates[:10]
    axes[2].barh([candidate.aptamer_id for candidate in top10], [candidate.composite_score for candidate in top10])
    axes[2].invert_yaxis()
    axes[2].set(title="Top candidates in rank order", xlabel="Enrichment ranking score")
    fig.tight_layout()
    try:
        fig.savefig(output_dir / "pipeline_summary.png", dpi=150)
    finally:
        plt.close(fig)


def export_results(ranked_candidates: list, target_features,
                   output_dir: Path, fmt: str = "csv") -> None:
    """Export ranked candidates to CSV and/or JSON."""
    if fmt not in ("csv", "json", "both"):
        raise ValueError("output.format must be csv, json, or both.")
    if not ranked_candidates:
        logger.warning("No candidates passed the filters. Writing empty results so previous winners are not mistaken for this run.")
        # Clear both result formats and the old plot, even if this run requested only CSV.
        fmt = "both"
        (output_dir / "pipeline_summary.png").unlink(missing_ok=True)

    records = [r.to_dict() for r in ranked_candidates]

    if fmt in ("csv", "both"):
        csv_path = output_dir / "ranked_candidates.csv"
        fieldnames = list(records[0].keys()) if records else [field.name for field in fields(RankedCandidate)]
        with open(csv_path, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(records)
        logger.info(f"Results exported to {csv_path}")

    if fmt in ("json", "both"):
        json_path = output_dir / "ranked_candidates.json"
        results = {
            "target": target_features.to_dict(),
            "n_candidates": len(ranked_candidates),
            "candidates": [r.to_dict() for r in ranked_candidates],
        }
        with open(json_path, "w") as f:
            json.dump(results, f, indent=2)
        logger.info(f"Results exported to {json_path}")


def run_pipeline(config: dict, stage: str = "all") -> dict:
    """Execute the aptamer identification pipeline.

    Args:
        config: Loaded pipeline configuration.
        stage: Which stage to run ("all" for full pipeline).

    Returns:
        Dictionary with pipeline results.
    """
    if stage not in STAGES:
        raise ValueError(f"Unknown pipeline stage: {stage}")
    results = {}
    output_config = config.get("output", {})
    output_dir = ensure_output_dir(output_config.get("directory", "output"))

    run_record = None
    if stage == "all":
        if output_config.get("format", "csv") not in ("csv", "json", "both"):
            raise ValueError("output.format must be csv, json, or both.")
        run_record = start_run_record(config, output_dir)

    needs_target = stage in ("target", "all")
    needs_library = stage in ("library", "structure", "scoring", "filtering", "all")
    needs_structure = stage == "structure"
    needs_scoring = stage in ("scoring", "filtering", "all")
    needs_filtering = stage in ("filtering", "all")

    # Station 4: Security Check validates the target for full pipeline runs.
    if needs_target:
        logger.info("=" * 60)
        logger.info("STAGE 1: Target Analysis")
        logger.info("=" * 60)
        target = analyze_target(config)
        results["target"] = target
        logger.info(f"Target: {target.name} ({target.length} residues)")

        if stage == "target":
            return results

    # Station 1 + 2: The Scanner + The Starting Line (real counts + CPM-ready candidates).
    if needs_library:
        logger.info("=" * 60)
        logger.info("STAGE 2: SELEX Count Ingestion")
        logger.info("=" * 60)
        candidates = generate_library(config)
        results["candidates"] = candidates
        logger.info(f"SELEX candidates: {len(candidates)} sequences retained")

        if stage == "library":
            return results

    # Structure prediction is now an optional annotation-only stage.
    if needs_structure:
        if "candidates" not in results:
            logger.error(
                "Library generation must run before structure prediction. "
                "Tip: run stage 'library' first."
            )
            return results

        logger.info("=" * 60)
        logger.info("STAGE 3: Structure Prediction")
        logger.info("=" * 60)
        structures = predict_structures(results["candidates"], config)
        results["structures"] = structures

        if stage == "structure":
            return results

    # Station 3: The Race Begins scoring from round trajectories.
    if needs_scoring:
        if "candidates" not in results:
            logger.error(
                "Library generation must run before scoring. "
                "Tip: run stage 'library' first."
            )
            return results

        logger.info("=" * 60)
        logger.info("STAGE 4: Enrichment Trajectory Scoring")
        logger.info("=" * 60)
        scores = score_binding(
            results["candidates"], config=config
        )
        results["binding_scores"] = scores

        if stage == "scoring":
            return results

    # Station 5: The Winning Bunch filtering + shortlist ranking.
    if needs_filtering:
        required = ["candidates", "binding_scores"]
        if not all(k in results for k in required):
            logger.error(
                "Previous stages must run before filtering. "
                "Tip: run through scoring first."
            )
            return results

        logger.info("=" * 60)
        logger.info("STAGE 5: Filtering & Ranking")
        logger.info("=" * 60)
        ranked = filter_and_rank(
            results["candidates"], results.get("structures", []),
            results["binding_scores"], config
        )
        results["ranked"] = ranked

    # Output
    if stage == "all" and "ranked" in results:
        logger.info("=" * 60)
        logger.info("STAGE 6: Results Export")
        logger.info("=" * 60)

        export_results(
            results["ranked"], results["target"], output_dir,
            fmt=output_config.get("format", "csv")
        )

        if output_config.get("generate_plots", True):
            generate_plots(results["ranked"], output_dir)

        logger.info("=" * 60)
        finish_run_record(
            run_record, output_dir, results["candidates"][0].round_order, len(results["ranked"])
        )
        logger.info("PIPELINE COMPLETE")
        logger.info(f"Top candidate: {results['ranked'][0].aptamer_id} "
                     f"(score={results['ranked'][0].composite_score:.4f})"
                     if results["ranked"] else "No candidates found")
        logger.info("=" * 60)

    return results


def main():
    parser = argparse.ArgumentParser(
        description="Aptamer Target Identification Pipeline",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python -m src.pipeline --config config/pipeline_config.yaml
  python -m src.pipeline --config config/pipeline_config.yaml --stage library
  python -m src.pipeline --config config/pipeline_config.yaml --top-n 100
        """
    )
    parser.add_argument(
        "--config", "-c", required=True,
        help="Path to pipeline YAML configuration file"
    )
    parser.add_argument(
        "--stage", "-s", choices=STAGES, default="all",
        help="Run a specific stage only (default: all)"
    )
    parser.add_argument(
        "--top-n", type=int, default=None,
        help="Override number of top candidates to return"
    )
    parser.add_argument(
        "--quiet", "-q", action="store_true",
        help="Suppress verbose output"
    )

    args = parser.parse_args()

    config = load_config(args.config)

    if args.top_n is not None:
        config.setdefault("filtering", {})["top_n"] = args.top_n

    verbose = not args.quiet and config.get("output", {}).get("verbose", True)
    setup_logging(verbose=verbose)

    logger.info("Aptamer Target Identification Pipeline")
    logger.info(f"Config: {args.config}")
    logger.info(f"Stage: {args.stage}")

    results = run_pipeline(config, stage=args.stage)

    if "ranked" in results and results["ranked"]:
        print(f"\nTop 5 candidates:")
        print(f"{'Rank':<6}{'ID':<14}{'Score':<10}{'log2E':<10}{'Guard':<10}{'Length':<8}{'GC':<8}")
        print("-" * 66)
        for r in results["ranked"][:5]:
            print(f"{r.rank:<6}{r.aptamer_id:<14}{r.composite_score:<10.4f}"
                  f"{r.log2_enrichment:<10.3f}{r.terminal_guardrail:<10.3f}"
                  f"{r.length:<8}{r.gc_content:<8.3f}")


if __name__ == "__main__":
    main()
