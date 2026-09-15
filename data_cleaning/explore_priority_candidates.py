from pathlib import Path

import pandas as pd


INPUT_FILE = Path("prioritised_hypothetical_with_go.csv")
OUTPUT_FILE = Path("top_10_priority_candidates.csv")
TOP_N = 10


def main() -> None:
    """Create a CSV containing the 10 unnamed proteins with the most computed GO terms."""

    if not INPUT_FILE.exists():
        raise FileNotFoundError(
            f"Could not find {INPUT_FILE}. Run parse_annotation_table.py first."
        )

    candidates = pd.read_csv(
        INPUT_FILE,
        dtype=str,
        keep_default_na=False
    )

    required_columns = [
        "Gene ID",
        "Product Description",
        "Computed GO Components",
        "Computed GO Functions",
        "Computed GO Processes",
        "EC numbers",
        "EC numbers from OrthoMCL",
        "GO component count",
        "GO function count",
        "GO process count",
        "Total computed GO count"
    ]

    missing_columns = [
        column for column in required_columns
        if column not in candidates.columns
    ]

    if missing_columns:
        raise ValueError(
            "The following required columns are missing:\n"
            + "\n".join(missing_columns)
        )

    count_columns = [
        "GO component count",
        "GO function count",
        "GO process count",
        "Total computed GO count"
    ]

    for column in count_columns:
        candidates[column] = pd.to_numeric(
            candidates[column],
            errors="coerce"
        ).fillna(0).astype(int)

    top_candidates = (
        candidates
        .sort_values(
            by=[
                "Total computed GO count",
                "GO function count",
                "GO process count",
                "Gene ID"
            ],
            ascending=[False, False, False, True]
        )
        .head(TOP_N)
        .copy()
    )

    output_columns = [
        "Gene ID",
        "Product Description",
        "GO component count",
        "GO function count",
        "GO process count",
        "Total computed GO count",
        "Computed GO Components",
        "Computed GO Functions",
        "Computed GO Processes",
        "EC numbers",
        "EC numbers from OrthoMCL"
    ]

    top_candidates[output_columns].to_csv(
        OUTPUT_FILE,
        index=False
    )

    print(f"Created: {OUTPUT_FILE.resolve()}")
    print(f"Number of candidates: {len(top_candidates)}")
    print("\nTop 10 candidates:")
    print(
        top_candidates[
            [
                "Gene ID",
                "GO component count",
                "GO function count",
                "GO process count",
                "Total computed GO count"
            ]
        ].to_string(index=False)
    )


if __name__ == "__main__":
    main()