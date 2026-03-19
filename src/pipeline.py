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
from pathlib import Path

from src.utils import load_config, setup_logging, ensure_output_dir
from src.sequence_generator import generate_library
from src.target_analyzer import analyze_target
from src.structure_predictor import predict_structures
from src.binding_scorer import score_binding
from src.filter_rank import filter_and_rank

logger = logging.getLogger("aptamer_pipeline")

STAGES = ["target", "library", "structure", "scoring", "filtering", "all"]


def generate_plots(ranked_candidates: list, output_dir: Path) -> None:
    """Generate summary visualization plots."""
    if not ranked_candidates:
        logger.warning("No candidates to plot.")
        return

    try:
        import pandas as pd
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        import seaborn as sns
    except Exception as exc:
        logger.warning(
            f"Plot generation skipped because plotting dependencies failed: {exc}. "
            "Tip: install matplotlib/seaborn or disable plotting in config."
        )
        return

    df = pd.DataFrame([r.to_dict() for r in ranked_candidates])

    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    fig.suptitle("Aptamer Pipeline — Top Candidates", fontsize=14, fontweight="bold")

    # 1. Composite score distribution
    sns.histplot(df["composite_score"], bins=20, ax=axes[0, 0], color="steelblue")
    axes[0, 0].set_title("Composite Score Distribution")
    axes[0, 0].set_xlabel("Composite Score")

    # 2. Binding vs Diversity scatter
    sns.scatterplot(data=df, x="diversity_score", y="binding_score",
                    hue="has_g_quadruplex", ax=axes[0, 1], alpha=0.7)
    axes[0, 1].set_title("Binding vs Sequence Diversity")
    axes[0, 1].set_xlabel("Diversity Score")
    axes[0, 1].set_ylabel("Binding Score")

    # 3. MFE distribution
    sns.histplot(df["mfe"], bins=20, ax=axes[1, 0], color="coral")
    axes[1, 0].set_title("Minimum Free Energy Distribution")
    axes[1, 0].set_xlabel("MFE (kcal/mol)")

    # 4. Score breakdown for top 10
    top10 = df.head(10)
    score_cols = ["binding_score", "diversity_score"]
    top10[score_cols].plot(kind="bar", stacked=True, ax=axes[1, 1],
                           color=["steelblue", "seagreen"])
    axes[1, 1].set_title("Score Breakdown — Top 10 Candidates")
    axes[1, 1].set_xlabel("Candidate")
    axes[1, 1].set_xticklabels(top10["aptamer_id"], rotation=45, ha="right")

    plt.tight_layout()
    plot_path = output_dir / "pipeline_summary.png"
    plt.savefig(plot_path, dpi=150)
    plt.close()
    logger.info(f"Plots saved to {plot_path}")


def export_results(ranked_candidates: list, target_features,
                   output_dir: Path, fmt: str = "csv") -> None:
    """Export ranked candidates to CSV and/or JSON."""
    if not ranked_candidates:
        logger.warning("No candidates to export.")
        return

    records = [r.to_dict() for r in ranked_candidates]

    if fmt in ("csv", "both"):
        csv_path = output_dir / "ranked_candidates.csv"
        fieldnames = list(records[0].keys())
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
    results = {}
    output_config = config.get("output", {})
    output_dir = ensure_output_dir(output_config.get("directory", "output"))

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
