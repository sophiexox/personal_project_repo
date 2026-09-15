from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# file paths
ROOT = Path(__file__).resolve().parents[1]

INPUT_FILE = ROOT / "Results" / "17_protnote_interpro2go_overlap_by_namespace.csv"
OUTPUT_FILE = ROOT / "Figures" / "18b_protnote_exact_agreement_by_namespace.png"

THRESHOLDS = [0.5, 0.7, 0.8, 0.9]

NAMESPACE_LABELS = {
    "MF": "Molecular Function",
    "CC": "Cellular Component",
    "BP": "Biological Process",
}

COLOURS = {
    "MF": "#6F58A8",
    "CC": "#5F9EA0",
    "BP": "#C36F85",
}


def main():
    if not INPUT_FILE.exists():
        raise FileNotFoundError(INPUT_FILE)

    df = pd.read_csv(INPUT_FILE)

    needed = {"Threshold", "Namespace", "Protein Exact Agreement Rate"}
    missing = needed - set(df.columns)

    if missing:
        raise KeyError(f"Missing columns: {', '.join(sorted(missing))}")

    df["Threshold"] = pd.to_numeric(df["Threshold"], errors="raise")
    df["Protein Exact Agreement Rate"] = pd.to_numeric(
        df["Protein Exact Agreement Rate"],
        errors="raise",
    )

    x = np.arange(len(THRESHOLDS))
    bar_width = 0.23

    fig, ax = plt.subplots(figsize=(11, 7))

    # plot one grouped bar series for each GO namespace
    for i, namespace in enumerate(NAMESPACE_LABELS):
        subset = (
            df[df["Namespace"] == namespace]
            .set_index("Threshold")
            .reindex(THRESHOLDS)
        )

        values = subset["Protein Exact Agreement Rate"].to_numpy() * 100
        positions = x + (i - 1) * bar_width

        bars = ax.bar(
            positions,
            values,
            width=bar_width,
            label=NAMESPACE_LABELS[namespace],
            color=COLOURS[namespace],
        )

        for bar, value in zip(bars, values):
            ax.text(
                bar.get_x() + bar.get_width() / 2,
                value + 1,
                f"{value:.1f}%",
                ha="center",
                va="bottom",
                fontsize=10,
            )

    ax.set_title(
        "ProtNote–InterPro2GO exact agreement by GO namespace",
        fontsize=19,
        pad=20,
    )
    ax.set_xlabel("ProtNote probability threshold", fontsize=13)
    ax.set_ylabel("Proteins with ≥1 exact GO match (%)", fontsize=13)
    ax.set_xticks(x)
    ax.set_xticklabels([f"≥{threshold:.1f}" for threshold in THRESHOLDS], fontsize=12)
    ax.set_ylim(0, 65)
    ax.tick_params(axis="y", labelsize=11)
    ax.grid(axis="y", alpha=0.18)
    ax.set_axisbelow(True)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.legend(frameon=False, fontsize=11, loc="upper right")

    plt.tight_layout()

    OUTPUT_FILE.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(OUTPUT_FILE, dpi=300, bbox_inches="tight")
    plt.show()


if __name__ == "__main__":
    main()


# GAI declaration
# GAI was used for lines 51-68 when arranging three namespace series as grouped
# bars around the same threshold positions. i asked how to offset several bar
# groups from a shared x position, and the example response used:
#
# x = np.arange(len(categories))
# positions = x + (group_index - 1) * bar_width
# ax.bar(positions, values, width=bar_width)
#
# i adapted that pattern for MF, CC and BP at each ProtNote threshold.

# i also used GAI for lines 70-78 to place percentage labels above each bar.
# the example calculated the centre with
# `bar.get_x() + bar.get_width() / 2` and passed that value to ax.text().
# i used the same approach for the namespace agreement percentages.