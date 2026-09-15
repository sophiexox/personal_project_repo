from pathlib import Path

import pandas as pd

# file paths
ROOT = Path(__file__).resolve().parents[1]

FUNGIDB_FILE = ROOT / "All_auris_genes_beta_fungidb_rel70.tsv"
MASTER_FILE = ROOT / "candida_master_mapping.xlsx"
RESULTS_DIR = ROOT / "Results"

EVIDENCE_COLUMNS = [
    "PFam ID",
    "CDD ID",
    "Gene3D ID",
    "Hamap ID",
    "NCBIfam ID",
    "PANTHER ID",
    "PRINTS ID",
    "SFLD ID",
    "PirSF ID",
    "Prositefamilies ID",
    "Smart ID",
    "Superfamily ID",
    "InterPro ID",
    "Computed GO Component IDs",
    "Computed GO Function IDs",
    "Computed GO Process IDs",
    "Curated GO Component IDs",
    "Curated GO Function IDs",
    "Curated GO Process IDs",
    "EC numbers",
    "EC numbers from OrthoMCL",
]


def has_value(series):
    values = series.fillna("").astype(str).str.strip()

    return ~(
        values.eq("")
        | values.str.upper().eq("N/A")
        | values.str.casefold().eq("not assigned")
    )


def clean_ids(series):
    return series.fillna("").astype(str).str.strip()


def main():
    for path in [FUNGIDB_FILE, MASTER_FILE]:
        if not path.exists():
            raise FileNotFoundError(path)

    fungidb = pd.read_csv(FUNGIDB_FILE, sep="\t", dtype=str)
    master = pd.read_excel(MASTER_FILE, sheet_name="Master", dtype=str)

    if "Gene ID" not in fungidb.columns:
        raise KeyError("Gene ID column not found in FungiDB data")

    if "VPDB Gene ID" not in master.columns:
        raise KeyError("VPDB Gene ID column not found in Master sheet")

    fungidb["Gene ID"] = clean_ids(fungidb["Gene ID"])
    master["VPDB Gene ID"] = clean_ids(master["VPDB Gene ID"])

    # find FungiDB records that never entered the master mapping
    master_ids = set(master.loc[master["VPDB Gene ID"].ne(""), "VPDB Gene ID"])
    unmatched = fungidb[~fungidb["Gene ID"].isin(master_ids)].copy()

    product = unmatched["Product Description"].fillna("").str.strip()
    unmatched["Is_Hypothetical"] = product.str.casefold().eq("hypothetical protein")

    # check conventional annotation evidence in the unmatched records
    found_evidence = [
        column for column in EVIDENCE_COLUMNS
        if column in unmatched.columns
    ]

    if not found_evidence:
        raise ValueError("No annotation evidence columns found")

    evidence_matrix = pd.concat(
        [has_value(unmatched[column]) for column in found_evidence],
        axis=1,
    )

    unmatched["Has_Conventional_Evidence"] = evidence_matrix.any(axis=1)
    unmatched["Annotation_Poor"] = (
        unmatched["Is_Hypothetical"]
        & ~unmatched["Has_Conventional_Evidence"]
    )

    RESULTS_DIR.mkdir(exist_ok=True)

    unmatched.to_csv(
        RESULTS_DIR / "09_all_fungidb_entries_missing_from_master.csv",
        index=False,
    )

    unmatched[unmatched["Annotation_Poor"]].to_csv(
        RESULTS_DIR / "09_annotation_poor_proteins_missing_from_master.csv",
        index=False,
    )

    hypothetical = int(unmatched["Is_Hypothetical"].sum())
    with_evidence = int(
        (
            unmatched["Is_Hypothetical"]
            & unmatched["Has_Conventional_Evidence"]
        ).sum()
    )
    annotation_poor = int(unmatched["Annotation_Poor"].sum())

    print(f"FungiDB entries absent from master: {len(unmatched)}")
    print(f"hypothetical among unmatched: {hypothetical}")
    print(f"hypothetical with evidence: {with_evidence}")
    print(f"annotation-poor among unmatched: {annotation_poor}")

    # a few extra checks that were useful when investigating the missing records
    if "Protein Length" in unmatched.columns:
        lengths = pd.to_numeric(unmatched["Protein Length"], errors="coerce")
        print(f"entries with protein length: {lengths.notna().sum()}")
        print(f"median protein length: {lengths.median():.1f}")

    print("\nevidence counts:")
    for column in found_evidence:
        print(f"{column}: {int(has_value(unmatched[column]).sum())}")


if __name__ == "__main__":
    main()


# GAI declaration
# GAI was used in this script for lines 69-70.
# i asked for a generated response on how to find rows in one dataframe whose IDs
# do not appear in another dataframe. i was given a generic example converting
# "your_master_dataframe['your_id_column']" to a set and then filtering
# "your_other_dataframe" with ~.isin(your_id_set). that was adapted to identify
# FungiDB Gene IDs that were absent from the VPDB Gene IDs in the Master sheet.

# GAI was also used in this script for lines 76-91.
# i asked how to apply the same missing-value check across whichever annotation
# columns were actually present in a dataframe, then create one boolean showing
# whether any evidence existed for each row. the response used a generic list
# comprehension over "your_evidence_columns", pd.concat(..., axis=1) and
# .any(axis=1). that was adapted to classify unmatched FungiDB records as having
# conventional evidence or being annotation-poor.

# GAI was also used in this script for lines 123-125.
# i asked for a generated response on how to convert a text column containing
# protein lengths into numbers while turning invalid values into missing values,
# then calculate a median. i was given a generic example using
# pd.to_numeric(your_dataframe["your_column"], errors="coerce") followed by
# .median(). that was adapted to inspect protein lengths among the unmatched
# FungiDB records.