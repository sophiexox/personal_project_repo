from pathlib import Path
import textwrap

import matplotlib.pyplot as plt
import pandas as pd

# file paths
ROOT = Path(__file__).resolve().parents[1]
RESULTS_DIR = ROOT / "Results"
FIGURES_DIR = ROOT / "Figures"

FLAGGED_FILE = RESULTS_DIR / "23_go2sum_flagged_outputs.csv"
SELECTED_FILE = RESULTS_DIR / "25_go2sum_selected_problematic_examples.csv"
FIGURE_PNG = FIGURES_DIR / "25_go2sum_problematic_examples.png"
FIGURE_PDF = FIGURES_DIR / "25_go2sum_problematic_examples.pdf"

# manually chosen examples from the flagged GO2SUM outputs
CATEGORIES = {
    "obsolete / placeholder": {
        "flag": "obsolete_placeholder",
        "explanation": "literal obsolete or malformed placeholder wording",
        "examples": [("B9J08_000672", "function"), ("B9J08_002843", "function")],
    },
    "photosynthesis / chloroplast": {
        "flag": "photosynthesis_or_chloroplast",
        "explanation": "plant or plastid language in this fungal protein set",
        "examples": [("B9J08_000070", "pathway"), ("B9J08_004645", "function")],
    },
    "metazoan / clinical": {
        "flag": "metazoan_or_clinical_context",
        "explanation": "human, animal or clinical wording needing taxonomic review",
        "examples": [("B9J08_005063", "subunit"), ("B9J08_000006", "function")],
    },
    "bacterial / viral": {
        "flag": "bacterial_or_viral_context",
        "explanation": "bacterial, viral or host-pathogen wording needing taxonomic review",
        "examples": [("B9J08_001949", "function"), ("B9J08_002031", "function")],
    },
    "extreme sequence-length wording": {
        "flag": "extreme_length_outlier",
        "explanation": "unusually long output by the per-output-type IQR rule",
        "examples": [("B9J08_001011", "function"), ("B9J08_003532", "function")],
    },
}


def shorten(text, width=190):
    # the CSV keeps the full description; only the figure copy is shortened
    return textwrap.shorten(str(text), width=width, placeholder=" …")


def main():
    if not FLAGGED_FILE.exists():
        raise FileNotFoundError(FLAGGED_FILE)

    flagged = pd.read_csv(FLAGGED_FILE)

    required = {
        "Protein",
        "Output_Type",
        "Description",
        "Word_Count",
        "Review_Reasons",
    }
    missing = required - set(flagged.columns)

    if missing:
        raise KeyError(f"Missing columns: {', '.join(sorted(missing))}")

    selected_rows = []

    for category, info in CATEGORIES.items():
        source_flag = info["flag"]

        # total number of outputs carrying this review flag
        available = flagged["Review_Reasons"].fillna("").str.contains(
            rf"(?:^|; ){source_flag}(?:;|$)",
            regex=True,
        ).sum()

        for protein, output_type in info["examples"]:
            match = flagged[
                flagged["Protein"].eq(protein)
                & flagged["Output_Type"].eq(output_type)
            ]

            if len(match) != 1:
                raise ValueError(
                    f"Expected one row for {protein} {output_type}, found {len(match)}"
                )

            row = match.iloc[0].copy()
            reasons = str(row["Review_Reasons"]).split("; ")

            if source_flag not in reasons:
                raise ValueError(
                    f"{protein} {output_type} does not carry {source_flag}"
                )

            row["Example_Category"] = category
            row["Category_Explanation"] = info["explanation"]
            row["Category_Source_Count"] = int(available)
            selected_rows.append(row)

    selected = pd.DataFrame(selected_rows)

    front = [
        "Example_Category",
        "Category_Explanation",
        "Category_Source_Count",
    ]
    selected = selected[
        front + [column for column in flagged.columns if column not in front]
    ]

    RESULTS_DIR.mkdir(exist_ok=True)
    FIGURES_DIR.mkdir(exist_ok=True)
    selected.to_csv(SELECTED_FILE, index=False)

    # build the compact table shown in the dissertation figure
    table_rows = []

    for row in selected.itertuples(index=False):
        table_rows.append([
            textwrap.fill(
                f"{row.Example_Category} (n={row.Category_Source_Count})\n"
                f"{row.Category_Explanation}",
                width=34,
            ),
            f"{row.Protein}\n{row.Output_Type}",
            textwrap.fill(shorten(row.Description), width=64),
            str(row.Word_Count),
        ])

    fig, ax = plt.subplots(figsize=(13, 10.2))
    ax.axis("off")

    table = ax.table(
        cellText=table_rows,
        colLabels=[
            "flag category",
            "protein / output",
            "original GO2SUM wording",
            "words",
        ],
        colWidths=[0.29, 0.14, 0.49, 0.08],
        cellLoc="left",
        colLoc="left",
        bbox=[0, 0.04, 1, 0.88],
    )
    table.auto_set_font_size(False)
    table.set_fontsize(8.5)
    table.scale(1, 2.2)

    for column in range(4):
        table[(0, column)].set_text_props(weight="bold")

    for row_number in range(1, len(table_rows) + 1):
        table[(row_number, 0)].set_fontsize(7.5)
        table[(row_number, 3)].get_text().set_ha("center")

    fig.suptitle(
        "Representative GO2SUM outputs flagged for manual review",
        fontsize=15,
        y=0.97,
    )
    fig.text(
        0.01,
        0.015,
        "Flags are lexical or length-based review prompts, not proof that a statement is false. "
        "Ellipses mark figure-only truncation; the CSV keeps the full original output.",
        fontsize=9,
    )

    fig.savefig(FIGURE_PNG, dpi=300, bbox_inches="tight")
    fig.savefig(FIGURE_PDF, bbox_inches="tight")
    plt.close(fig)

    print(f"selected examples: {len(selected)}")


if __name__ == "__main__":
    main()


# GAI declaration
# GAI was used for lines 68-99 when selecting manually chosen Protein/Output_Type
# combinations from the larger flagged-results table. i asked how to match on
# two dataframe columns at once and check that exactly one source row had been
# returned. the example response used:
#
# match = df[
#     df["your_id"].eq(identifier)
#     & df["your_type"].eq(output_type)
# ]
#
# if len(match) != 1:
#     raise ValueError(...)
#
# i adapted this to the selected GO2SUM examples and added a second check that
# each row still carried the review flag expected for its category.

# GAI was also used for lines 116-154 when constructing the matplotlib table.
# i asked how to display long dataframe text in a figure without allowing the
# descriptions to run outside the cells. the example used textwrap.fill() for
# wrapping, ax.table() with manually chosen column widths, and direct cell access
# such as table[(row, column)] to alter font size and alignment. i adapted this
# into the four-column figure while keeping the full unshortened descriptions in
# the accompanying CSV.