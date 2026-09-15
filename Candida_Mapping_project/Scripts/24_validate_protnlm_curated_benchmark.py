from pathlib import Path
import re

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# file paths
ROOT = Path(__file__).resolve().parents[1]
PROTNLM_DIR = ROOT / "Data" / "raw" / "protnlm"
RESULTS_DIR = ROOT / "Results"
FIGURES_DIR = ROOT / "Figures"

BENCHMARK_FILE = PROTNLM_DIR / "ProtNLM_curated_benchmark_94_predictions.csv"
ANNOTATION_POOR_FILE = PROTNLM_DIR / "ProtNLM_annotation_poor_521_predictions.csv"

REVIEW_FILE = RESULTS_DIR / "24_protnlm_curated_benchmark_manual_review.csv"
SUMMARY_FILE = RESULTS_DIR / "24_protnlm_benchmark_score_summary.csv"
COMPARISON_PNG = FIGURES_DIR / "24_protnlm_benchmark_comparison.png"
COMPARISON_PDF = FIGURES_DIR / "24_protnlm_benchmark_comparison.pdf"

CLASSIFIED_REVIEW_FILE = RESULTS_DIR / "24_protnlm_curated_benchmark_manual_review_first_pass.csv"
ASSESSMENT_SUMMARY_FILE = RESULTS_DIR / "24_protnlm_score_by_assessment_summary.csv"
ASSESSMENT_PNG = FIGURES_DIR / "24_protnlm_score_by_assessment.png"
ASSESSMENT_PDF = FIGURES_DIR / "24_protnlm_score_by_assessment.pdf"

ASSESSMENT_ORDER = ["specific", "related", "discordant", "uncertain"]


def parse_header(header):
    # pull a compact protein ID, gene name and existing annotation from the FASTA header
    text = str(header)

    protein_match = re.search(r"B9J08_\d+", text)
    protein_id = protein_match.group(0) if protein_match else ""

    gene_match = re.search(r"B9J08_\d+\s+([^\s]+)", text)
    gene_name = gene_match.group(1) if gene_match else ""

    marker = "Verified ORF;"
    annotation = text.split(marker, 1)[1].strip() if marker in text else text

    return protein_id, gene_name, annotation


def score_summary(name, scores):
    # one compact row describing a rank-1 score distribution
    return {
        "dataset": name,
        "n_proteins": len(scores),
        "mean_score": scores.mean(),
        "median_score": scores.median(),
        "minimum_score": scores.min(),
        "maximum_score": scores.max(),
        "score_ge_0.2": int((scores >= 0.2).sum()),
        "score_ge_0.5": int((scores >= 0.5).sum()),
        "score_ge_0.8": int((scores >= 0.8).sum()),
        "score_ge_0.9": int((scores >= 0.9).sum()),
    }


def make_comparison_figure(benchmark_scores, annotation_poor_scores):
    bins = np.arange(0, 1.05, 0.05)

    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    # percentages make the 94- and 521-protein groups directly comparable
    axes[0].hist(
        annotation_poor_scores,
        bins=bins,
        weights=np.full(len(annotation_poor_scores), 100 / len(annotation_poor_scores)),
        alpha=0.55,
        label=f"annotation-poor proteins (n={len(annotation_poor_scores)})",
    )
    axes[0].hist(
        benchmark_scores,
        bins=bins,
        weights=np.full(len(benchmark_scores), 100 / len(benchmark_scores)),
        alpha=0.55,
        label=f"curated benchmark (n={len(benchmark_scores)})",
    )
    axes[0].axvline(0.2, linestyle="--", linewidth=1.2, label="score = 0.2")
    axes[0].set_xlabel("ProtNLM top-ranked prediction score")
    axes[0].set_ylabel("proteins (%)")
    axes[0].set_title("a  Normalised score distributions", loc="left")
    axes[0].legend()

    axes[1].boxplot(
        [benchmark_scores, annotation_poor_scores],
        tick_labels=["curated benchmark", "annotation-poor"],
        showfliers=True,
    )
    axes[1].set_ylabel("ProtNLM top-ranked prediction score")
    axes[1].set_title("b  Score comparison between protein groups", loc="left")

    fig.tight_layout()
    fig.savefig(COMPARISON_PNG, dpi=300, bbox_inches="tight")
    fig.savefig(COMPARISON_PDF, bbox_inches="tight")
    plt.close(fig)


def make_assessment_figure(classified):
    # keep each manual category separate and show every reviewed protein
    scores = [
        classified.loc[
            classified["manual_assessment"] == category,
            "protnlm_score",
        ].dropna()
        for category in ASSESSMENT_ORDER
    ]

    fig, ax = plt.subplots(figsize=(8, 5))
    ax.boxplot(scores, tick_labels=ASSESSMENT_ORDER, showfliers=False)

    rng = np.random.default_rng(24)

    for position, values in enumerate(scores, start=1):
        jitter = rng.normal(position, 0.055, len(values))
        ax.scatter(jitter, values, alpha=0.7, s=28, zorder=3)

    ax.set_xlabel("first-pass assessment category")
    ax.set_ylabel("ProtNLM top-ranked prediction score")
    ax.set_title("ProtNLM score by first-pass assessment")

    fig.tight_layout()
    fig.savefig(ASSESSMENT_PNG, dpi=300, bbox_inches="tight")
    fig.savefig(ASSESSMENT_PDF, bbox_inches="tight")
    plt.close(fig)


def main():
    for path in [BENCHMARK_FILE, ANNOTATION_POOR_FILE]:
        if not path.exists():
            raise FileNotFoundError(path)

    RESULTS_DIR.mkdir(exist_ok=True)
    FIGURES_DIR.mkdir(exist_ok=True)

    benchmark = pd.read_csv(BENCHMARK_FILE)
    annotation_poor = pd.read_csv(ANNOTATION_POOR_FILE)

    needed = {"Protein", "Rank", "Prediction", "Score"}
    for name, frame in [("benchmark", benchmark), ("annotation-poor", annotation_poor)]:
        missing = needed - set(frame.columns)
        if missing:
            raise KeyError(f"Missing {name} columns: {', '.join(sorted(missing))}")

    benchmark["Rank"] = pd.to_numeric(benchmark["Rank"], errors="raise")
    benchmark["Score"] = pd.to_numeric(benchmark["Score"], errors="coerce")
    annotation_poor["Rank"] = pd.to_numeric(annotation_poor["Rank"], errors="raise")
    annotation_poor["Score"] = pd.to_numeric(annotation_poor["Score"], errors="coerce")

    # benchmark should contain 94 proteins with 10 predictions each
    if benchmark["Protein"].nunique() != 94:
        raise ValueError("Expected 94 curated benchmark proteins")

    if (benchmark.groupby("Protein").size() != 10).any():
        raise ValueError("Every curated benchmark protein should have 10 predictions")

    benchmark_top = benchmark[benchmark["Rank"] == 1].copy()
    annotation_poor_top = annotation_poor[annotation_poor["Rank"] == 1].copy()

    # create a score-blind table for manual biological review
    parsed = benchmark_top["Protein"].apply(parse_header)

    review = pd.DataFrame({
        "protein_id": parsed.str[0],
        "gene_name": parsed.str[1],
        "existing_annotation": parsed.str[2],
        "protnlm_top_prediction": benchmark_top["Prediction"].to_numpy(),
        "protnlm_score": benchmark_top["Score"].to_numpy(),
        "manual_assessment": "",
        "review_notes": "",
    })
    review.to_csv(REVIEW_FILE, index=False)

    benchmark_scores = benchmark_top["Score"].dropna()
    annotation_poor_scores = annotation_poor_top["Score"].dropna()

    summary = pd.DataFrame([
        score_summary("curated benchmark", benchmark_scores),
        score_summary("annotation-poor", annotation_poor_scores),
    ])
    summary.to_csv(SUMMARY_FILE, index=False)

    make_comparison_figure(benchmark_scores, annotation_poor_scores)

    print(f"curated benchmark proteins: {len(benchmark_scores)}")
    print(f"annotation-poor proteins: {len(annotation_poor_scores)}")
    print(f"benchmark median score: {benchmark_scores.median():.4f}")
    print(f"annotation-poor median score: {annotation_poor_scores.median():.4f}")

    # the manually completed first-pass file is analysed separately
    if not CLASSIFIED_REVIEW_FILE.exists():
        print("first-pass review file not found; assessment analysis skipped")
        return

    classified = pd.read_csv(CLASSIFIED_REVIEW_FILE)
    classified["manual_assessment"] = (
        classified["manual_assessment"].fillna("").str.strip().str.casefold()
    )
    classified["protnlm_score"] = pd.to_numeric(
        classified["protnlm_score"],
        errors="coerce",
    )

    found = set(classified["manual_assessment"]) - {""}
    unexpected = found - set(ASSESSMENT_ORDER)

    if unexpected:
        raise ValueError(f"Unexpected assessment categories: {sorted(unexpected)}")

    assessment_summary = (
        classified[classified["manual_assessment"].ne("")]
        .groupby("manual_assessment")["protnlm_score"]
        .agg(n="count", mean="mean", median="median", min="min", max="max")
        .reindex(ASSESSMENT_ORDER)
        .reset_index()
    )
    assessment_summary.to_csv(ASSESSMENT_SUMMARY_FILE, index=False)

    make_assessment_figure(classified)

    print("\nassessment score summary:")
    print(assessment_summary.to_string(index=False))


if __name__ == "__main__":
    main()


# GAI declaration
# GAI was used for lines 30-43 when extracting information from the long CGD
# FASTA-style protein headers. i asked how a regular expression could find an
# identifier matching B9J08 followed by digits and then capture the next word as
# a gene name. the example response used code similar to:
#
# id_match = re.search(r"YOUR_PREFIX_\d+", header)
# name_match = re.search(r"YOUR_PREFIX_\d+\s+([^\s]+)", header)
#
# identifier = id_match.group(0) if id_match else ""
# name = name_match.group(1) if name_match else ""
#
# i adapted this to the B8441 identifiers. GAI also suggested using
# text.split("your marker", 1)[1] when i asked how to retain everything after a
# known phrase in the header, which is used here for the existing annotation.

# lines 62-99 were written with GAI support when comparing score distributions
# from groups with very different sample sizes. i asked how histograms for 94 and
# 521 observations could be expressed as percentages rather than raw counts.
# ChatGPT suggested assigning every observation a weight of 100 divided by the
# number of observations:
#
# weights = np.full(len(scores), 100 / len(scores))
# ax.hist(scores, bins=bins, weights=weights)
#
# i adapted this for both datasets and combined the percentage histograms with a
# boxplot so the overall distributions could be compared without the larger
# annotation-poor group dominating simply because it contained more proteins.

# GAI was also used for lines 102-128 when visualising ProtNLM scores after the
# manual first-pass assessment. i wanted to retain a boxplot for each assessment
# group while also showing every individual protein. the example response used a
# small random horizontal offset around each category position:
#
# rng = np.random.default_rng(1)
# jitter = rng.normal(position, 0.05, len(values))
# ax.scatter(jitter, values)
#
# i adapted this with a fixed seed so the point positions remain reproducible.
# the manual categories themselves were assigned separately from the model score;
# this plot is only used afterwards to examine how the score distributions relate
# to those biological review categories.

# finally, lines 213-220 use a pandas aggregation pattern suggested by GAI. i
# asked how to calculate the number, mean, median, minimum and maximum score for
# each manual assessment category in one operation. the example used:
#
# summary = (
#     df.groupby("your_category")["your_score"]
#     .agg(n="count", mean="mean", median="median", min="min", max="max")
#     .reset_index()
# )
#
# i adapted this to manual_assessment and protnlm_score, then reindexed the table
# so the categories always appear in the intended specific/related/discordant/
# uncertain order.