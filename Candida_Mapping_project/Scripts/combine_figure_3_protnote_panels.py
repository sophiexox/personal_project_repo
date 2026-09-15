from pathlib import Path

import matplotlib.pyplot as plt

# file paths
ROOT = Path(__file__).resolve().parents[1]
FIGURES_DIR = ROOT / "Figures"

PANEL_FILES = [
    FIGURES_DIR / "18a_protnote_exact_agreement_by_threshold.png",
    FIGURES_DIR / "18b_protnote_exact_agreement_by_namespace.png",
]

OUTPUT_PNG = FIGURES_DIR / "Figure_3_ProtNote_agreement.png"
OUTPUT_PDF = FIGURES_DIR / "Figure_3_ProtNote_agreement.pdf"


def main():
    for path in PANEL_FILES:
        if not path.exists():
            raise FileNotFoundError(path)

    images = [plt.imread(path) for path in PANEL_FILES]

    fig, axes = plt.subplots(1, 2, figsize=(12.8, 5.4))

    for ax, image, label in zip(axes, images, "AB"):
        ax.imshow(image)
        ax.axis("off")
        ax.text(
            0.02,
            0.98,
            label,
            transform=ax.transAxes,
            ha="left",
            va="top",
            fontsize=18,
            fontweight="bold",
            bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.9},
        )

    fig.tight_layout(pad=0.5)
    FIGURES_DIR.mkdir(exist_ok=True)

    fig.savefig(OUTPUT_PNG, dpi=600, bbox_inches="tight", facecolor="white")
    fig.savefig(OUTPUT_PDF, bbox_inches="tight", facecolor="white")
    plt.close(fig)


if __name__ == "__main__":
    main()


# GAI declaration
# GAI was used for lines 27-40 when placing two existing PNG plots into one
# matplotlib figure and adding A/B panel labels. i adapted an example using
# zip(), imshow() and ax.text() to assemble the final dissertation figure.