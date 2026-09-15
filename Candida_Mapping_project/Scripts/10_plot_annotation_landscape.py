from pathlib import Path

import matplotlib.pyplot as plt

# file paths
ROOT = Path(__file__).resolve().parents[1]
FIGURE_DIR = ROOT / "Figures"

# final counts from the release 70 analysis
TOTAL_PROTEINS = 5417
HYPOTHETICAL = 860
WITH_EVIDENCE = 334
ANNOTATION_POOR = 526


def draw_stacked_bar(ax, labels, counts, title):
    # draw one horizontal stacked bar and label each section
    left = 0
    total = sum(counts)

    for label, count in zip(labels, counts):
        ax.barh([""], [count], left=left)

        ax.text(
            left + count / 2,
            0,
            f"{label}\n{count:,}\n({count / total * 100:.1f}%)",
            ha="center",
            va="center",
            fontsize=11,
        )

        left += count

    ax.set_title(title)
    ax.set_xlabel("Number of proteins")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)


def main():
    FIGURE_DIR.mkdir(exist_ok=True)

    named = TOTAL_PROTEINS - HYPOTHETICAL

    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    draw_stacked_bar(
        axes[0],
        ["Informatively named", "Hypothetical"],
        [named, HYPOTHETICAL],
        "A. Current protein annotation",
    )

    draw_stacked_bar(
        axes[1],
        ["Conventional evidence", "Annotation-poor"],
        [WITH_EVIDENCE, ANNOTATION_POOR],
        "B. Evidence among hypothetical proteins",
    )

    fig.suptitle("Annotation landscape of Candida auris B8441", fontsize=15)
    fig.text(0.5, 0.92, "Beta FungiDB release 70", ha="center", fontsize=10)

    plt.tight_layout(rect=[0, 0, 1, 0.88])

    # save both versions used for the dissertation and repository
    plt.savefig(
        FIGURE_DIR / "fungidb_rel70_annotation_landscape.png",
        dpi=300,
        bbox_inches="tight",
    )
    plt.savefig(
        FIGURE_DIR / "fungidb_rel70_annotation_landscape.pdf",
        bbox_inches="tight",
    )

    plt.show()


if __name__ == "__main__":
    main()


# GAI declaration
# GAI was used for lines 16-38 when writing the repeated plotting step as one
# function. i asked how two similar stacked horizontal bar plots could be drawn
# without duplicating the same matplotlib code twice. the example response used
# a function taking an axis, a list of labels and a list of counts, with a
# running "left" value passed into ax.barh(). i adapted this into
# draw_stacked_bar() for the two annotation panels.

# lines 24-30 also came from a separate GAI query about label placement. i asked
# how to put the count and percentage in the middle of each section of a stacked
# horizontal bar. ChatGPT suggested calculating the x position as
# "left + value / 2" and the percentage as "value / total * 100". i used those
# calculations for the text placed inside each bar segment.