from pathlib import Path

import pandas as pd

# file path
ROOT = Path(__file__).resolve().parents[1]
FUNGIDB_FILE = ROOT / "All_auris_genes_beta_fungidb_rel70.tsv"


def missing_mask(series):
    # treat blank, N/A and not assigned as missing annotation
    values = (
        series.fillna("")
        .astype(str)
        .str.strip()
    )

    return (
        values.eq("")
        | values.str.upper().eq("N/A")
        | values.str.casefold().eq("not assigned")
    )


def main():
    if not FUNGIDB_FILE.exists():
        raise FileNotFoundError(FUNGIDB_FILE)

    df = pd.read_csv(FUNGIDB_FILE, sep="\t", dtype=str)

    if "Product Description" not in df.columns:
        raise KeyError("Product Description column not found")

    product = df["Product Description"].fillna("").str.strip()
    hypothetical = product.str.casefold().eq("hypothetical protein")

    print(f"genes: {len(df)}")
    print(f"hypothetical proteins: {hypothetical.sum()}")
    print(f"non-hypothetical proteins: {(~hypothetical).sum()}")

    # find annotation columns without hard-coding every column name
    keywords = [
        "Pfam",
        "InterPro",
        "PANTHER",
        "Superfamily",
        "SMART",
        "GO",
        "Curated",
        "Computed",
        "Ortholog",
        "Paralog",
    ]

    annotation_cols = [
        column
        for column in df.columns
        if any(word.casefold() in column.casefold() for word in keywords)
    ]

    print("\nannotation columns:")
    for column in annotation_cols:
        present = (~missing_mask(df[column])).sum()
        print(f"{column}: {present} present")

    # check how much evidence already exists for hypothetical proteins
    evidence_words = ["Pfam", "InterPro", "PANTHER", "Superfamily", "SMART", "Computed GO"]
    evidence_cols = [
        column
        for column in df.columns
        if any(word.casefold() in column.casefold() for word in evidence_words)
    ]

    hypothetical_df = df[hypothetical]

    print("\nhypothetical protein evidence:")
    for column in evidence_cols:
        count = (~missing_mask(hypothetical_df[column])).sum()
        print(f"{column}: {count}")


if __name__ == "__main__":
    main()


# GAI declaration
# GAI was used in this script for lines 10-22.
# i asked for a generated response on how to avoid repeating the same checks for
# blank values, "N/A" and "not assigned" every time i inspected an annotation
# column. i was given a generic example which created a helper function using
# "your_series".fillna("").str.strip() and returned a boolean mask combining the
# different missing-value conditions. that was adapted into missing_mask() here.

# GAI was also used in this script for lines 55-72.
# i asked for a generated response on how to find dataframe columns containing
# any word from a list of annotation keywords without manually specifying every
# column name. i was given a generic list comprehension using
# "your_dataframe".columns and any(keyword in column for keyword in
# "your_keywords"). that approach was adapted to identify the annotation_cols
# and evidence_cols used in this script.