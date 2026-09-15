from pathlib import Path

import pandas as pd

# file paths
ROOT = Path(__file__).resolve().parents[1]
FUNGIDB_FILE = ROOT / "All_auris_genes_beta_fungidb_rel70.tsv"
RESULTS_DIR = ROOT / "Results"

# annotation columns used as conventional evidence
EVIDENCE_SOURCES = {
    "PFam": "PFam ID",
    "CDD": "CDD ID",
    "Gene3D": "Gene3D ID",
    "Hamap": "Hamap ID",
    "NCBIfam": "NCBIfam ID",
    "PANTHER": "PANTHER ID",
    "PRINTS": "PRINTS ID",
    "SFLD": "SFLD ID",
    "PirSF": "PirSF ID",
    "Prosite": "Prositefamilies ID",
    "SMART": "Smart ID",
    "Superfamily": "Superfamily ID",
    "InterPro": "InterPro ID",
}

GO_COLUMNS = [
    "Computed GO Component IDs",
    "Computed GO Function IDs",
    "Computed GO Process IDs",
]


def has_value(series):
    # True where an annotation field contains useful information
    values = series.fillna("").astype(str).str.strip()

    return ~(
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

    # select proteins still labelled as hypothetical
    product = df["Product Description"].fillna("").str.strip().str.casefold()
    hypothetical = df[product.eq("hypothetical protein")].copy()

    # add one boolean evidence flag for each annotation source
    evidence_flags = []

    for name, column in EVIDENCE_SOURCES.items():
        if column not in hypothetical.columns:
            raise KeyError(f"Missing evidence column: {column}")

        flag = f"Has_{name}"
        hypothetical[flag] = has_value(hypothetical[column])
        evidence_flags.append(flag)

    # combine the three computed GO fields into one GO evidence flag
    for column in GO_COLUMNS:
        if column not in hypothetical.columns:
            raise KeyError(f"Missing GO column: {column}")

    hypothetical["Has_Any_GO"] = hypothetical[GO_COLUMNS].apply(
        lambda column: has_value(column)
    ).any(axis=1)
    evidence_flags.append("Has_Any_GO")

    # proteins with no selected evidence are the annotation-poor set
    hypothetical["Has_Any_Conventional_Evidence"] = (
        hypothetical[evidence_flags].any(axis=1)
    )
    hypothetical["Annotation_Poor"] = ~hypothetical["Has_Any_Conventional_Evidence"]

    annotation_poor = hypothetical[hypothetical["Annotation_Poor"]].copy()
    with_evidence = hypothetical[
        hypothetical["Has_Any_Conventional_Evidence"]
    ].copy()

    # count the most common combinations of evidence types
    evidence_patterns = (
        hypothetical[evidence_flags]
        .value_counts()
        .reset_index(name="Protein_Count")
    )

    # smaller table for the proteins carried forward to AI analysis
    candidate_columns = [
        "Gene ID",
        "source_id",
        "Product Description",
        "Protein Length",
        "Ortholog count",
        "Paralog count",
        *EVIDENCE_SOURCES.values(),
        *GO_COLUMNS,
    ]
    candidate_columns = [
        column
        for column in candidate_columns
        if column in annotation_poor.columns
    ]
    candidates = annotation_poor[candidate_columns].copy()

    RESULTS_DIR.mkdir(exist_ok=True)

    outputs = {
        "fungidb_rel70_all_hypothetical_proteins.csv": hypothetical,
        "fungidb_rel70_annotation_poor_hypothetical_proteins.csv": annotation_poor,
        "fungidb_rel70_hypothetical_with_evidence.csv": with_evidence,
        "fungidb_rel70_hypothetical_evidence_patterns.csv": evidence_patterns,
        "fungidb_rel70_annotation_poor_candidate_list.csv": candidates,
    }

    for filename, data in outputs.items():
        data.to_csv(RESULTS_DIR / filename, index=False)

    total = len(hypothetical)
    evidence_count = int(
        hypothetical["Has_Any_Conventional_Evidence"].sum()
    )
    poor_count = len(annotation_poor)

    print(f"hypothetical proteins: {total}")
    print(f"with conventional evidence: {evidence_count}")
    print(f"annotation-poor: {poor_count}")
    print(f"annotation-poor percentage: {poor_count / total * 100:.2f}%")

    print("\nevidence counts:")
    for flag in evidence_flags:
        count = int(hypothetical[flag].sum())
        print(f"{flag}: {count} ({count / total * 100:.2f}%)")


if __name__ == "__main__":
    main()


# GAI declaration
# GAI was used in this script for lines 58-68.
# i asked for a generated response on how to create several boolean columns from
# a dictionary that maps short evidence names to dataframe column names, instead
# of writing the same assignment repeatedly. i was given a generic example using
# "your_mapping".items(), creating a new column such as f"Has_{name}", and
# appending the new column names to a list. that approach was adapted to create
# the PFam, CDD, Gene3D, InterPro and other evidence flags in this script.

# GAI was also used in this script for lines 70-78.
# i asked how to check several GO annotation columns with the same missing-value
# function and then return one True/False result per protein if any of the
# columns contained a value. the response used a generic example with
# "your_dataframe[your_columns].apply(your_function).any(axis=1)". that was
# adapted to the three computed GO columns and used to create Has_Any_GO.

# GAI was also used in this script for lines 80-90.
# i asked for a generated response on how to classify rows as having any evidence
# across several boolean columns, then create separate dataframes for rows with
# and without evidence. i was given a generic pandas example using
# "your_dataframe[your_boolean_columns].any(axis=1)" followed by boolean
# filtering. that was adapted to define Has_Any_Conventional_Evidence and the
# annotation-poor protein subset.

# GAI was also used in this script for lines 92-96.
# i asked how to count the most common combinations of True/False values across
# several evidence columns. the generated example used
# "your_dataframe[your_columns].value_counts().reset_index(name='Count')".
# that was adapted to create the evidence_patterns table used in the results.