from pathlib import Path

import matplotlib.pyplot as plt

# file paths
ROOT = Path(__file__).resolve().parents[1]
FIGURES_DIR = ROOT / "Figures"

PANEL_FILES = [
    FIGURES_DIR / "22_protnlm_score_vs_protein_length.png",
    FIGURES_DIR / "22_protnlm_score_vs_ortholog_count.png",
    FIGURES_DIR / "22_protnlm_score_vs_lexical_beam_coherence.png",
]

OUTPUT_PNG = FIGURES_DIR / "Figure_2_ProtNLM_diagnostics.png"
OUTPUT_PDF = FIGURES_DIR / "Figure_2_ProtNLM_diagnostics.pdf"


def main():
    for path in PANEL_FILES:
        if not path.exists():
            raise FileNotFoundError(path)

    images = [plt.imread(path) for path in PANEL_FILES]

    fig, axes = plt.subplots(1, 3, figsize=(14.5, 4.8))

    for ax, image, label in zip(axes, images, "ABC"):
        ax.imshow(image)
        ax.axis("off")
        ax.text(
            0.02,
            0.98,
            label,
            transform=ax.transAxes,
            ha="left",
            va="top",
            fontsize=16,
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
# GAI was used for lines 24-41 when combining three existing PNG plots into one
# matplotlib figure and adding A/B/C panel labels. i asked how to loop through
# several images and axes together, and adapted an example using zip(), imshow()
# and ax.text() for the final dissertation figure.