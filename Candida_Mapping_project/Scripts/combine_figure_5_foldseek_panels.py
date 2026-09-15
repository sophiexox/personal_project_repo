from pathlib import Path

import matplotlib.pyplot as plt

# file paths
ROOT = Path(__file__).resolve().parents[1]
FIGURES_DIR = ROOT / "Figures"

PANEL_FILES = [
    FIGURES_DIR / "B9J08_004442_Foldseek_Seipin_alignment.png",
    FIGURES_DIR / "B9J08_003070_foldseek.png",
]

OUTPUT_PNG = FIGURES_DIR / "Figure_5_Foldseek_structural_comparisons.png"
OUTPUT_PDF = FIGURES_DIR / "Figure_5_Foldseek_structural_comparisons.pdf"


def main():
    for path in PANEL_FILES:
        if not path.exists():
            raise FileNotFoundError(path)

    images = [plt.imread(path) for path in PANEL_FILES]

    fig, axes = plt.subplots(1, 2, figsize=(12.8, 5.8), facecolor="black")

    for ax, image, label in zip(axes, images, "AB"):
        ax.set_facecolor("black")
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
            bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.95},
        )

    fig.tight_layout(pad=0.3)
    FIGURES_DIR.mkdir(exist_ok=True)

    fig.savefig(OUTPUT_PNG, dpi=600, bbox_inches="tight", facecolor="black")
    fig.savefig(OUTPUT_PDF, bbox_inches="tight", facecolor="black")
    plt.close(fig)


if __name__ == "__main__":
    main()


# GAI declaration
# GAI was used for lines 27-40 when placing two existing structural PNG renders
# side by side on a black matplotlib figure and adding A/B panel labels. i
# adapted an example using zip(), imshow() and ax.text() to assemble the final
# dissertation figure without recreating the Foldseek results.