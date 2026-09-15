from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd

# file paths
ROOT = Path(__file__).resolve().parents[1]

INPUT_FILE = ROOT / "Results" / "16_protnote_interpro2go_exact_overlap_summary.csv"
FIGURE_DIR = ROOT / "Figures"

PNG_OUTPUT = FIGURE_DIR / "18a_protnote_exact_agreement_by_threshold.png"
PDF_OUTPUT = FIGURE_DIR / "18a_protnote_exact_agreement_by_threshold.pdf"


def main():
    if not INPUT_FILE.exists():
        raise FileNotFoundError(INPUT_FILE)

    df = pd.read_csv(INPUT_FILE).sort_values("Threshold").copy()

    df["Agreement_Percent"] = df["Protein Exact Agreement Rate"] * 100
    df["Coverage_Percent"] = (
        df["Proteins With ProtNote Predictions"]
        / df["Benchmark Proteins"]
        * 100
    )

    fig, ax = plt.subplots(figsize=(8.5, 5.5))

    # exact agreement and prediction coverage across thresholds
    ax.plot(
        df["Threshold"],
        df["Agreement_Percent"],
        marker="o",
        linewidth=2.8,
        markersize=8,
        color="#8E5EA2",
        label="Proteins with ≥1 exact GO match",
    )
    ax.plot(
        df["Threshold"],
        df["Coverage_Percent"],
        marker="s",
        linewidth=2.2,
        markersize=7,
        linestyle="--",
        color="#C06C84",
        label="Proteins retaining ≥1 ProtNote prediction",
    )

    # label each point with its percentage
    for row in df.itertuples(index=False):
        ax.text(
            row.Threshold,
            row.Agreement_Percent + 2.2,
            f"{row.Agreement_Percent:.1f}%",
            ha="center",
            va="bottom",
            fontsize=10,
        )
        ax.text(
            row.Threshold,
            row.Coverage_Percent - 4.0,
            f"{row.Coverage_Percent:.1f}%",
            ha="center",
            va="top",
            fontsize=9,
        )

    ax.set_title(
        "ProtNote–InterPro2GO exact agreement across probability thresholds",
        fontsize=14,
        pad=14,
    )
    ax.set_xlabel("ProtNote probability threshold", fontsize=11)
    ax.set_ylabel("Proteins (%)", fontsize=11)
    ax.set_xticks(df["Threshold"])
    ax.set_xticklabels([f"≥{value:.1f}" for value in df["Threshold"]])
    ax.set_ylim(0, 105)
    ax.grid(axis="y", alpha=0.2)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.legend(frameon=False, loc="lower left")

    plt.tight_layout()

    FIGURE_DIR.mkdir(exist_ok=True)
    plt.savefig(PNG_OUTPUT, dpi=300, bbox_inches="tight")
    plt.savefig(PDF_OUTPUT, bbox_inches="tight")
    plt.show()


if __name__ == "__main__":
    main()


# GAI declaration
# GAI was used for lines 22-27. i asked how to convert a proportion column to a
# percentage and calculate prediction coverage as one count divided by another.
# the example response used:
#
# df["percent"] = df["rate"] * 100
# df["coverage"] = df["predicted"] / df["total"] * 100
#
# i adapted those calculations to the exact agreement and benchmark coverage
# columns produced by script 16.

# i also used GAI for lines 53-69 when adding percentage labels to both plotted
# series. i asked how to loop through dataframe rows and position text slightly
# above one line and below another. the example used ax.text(x, y + offset, ...)
# inside a row loop, which i adapted for Agreement_Percent and Coverage_Percent.