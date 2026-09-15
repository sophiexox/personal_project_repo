from pathlib import Path
from datetime import datetime
import shutil

import pandas as pd

# file paths
ROOT = Path(__file__).resolve().parents[1]

FUNGIDB_FILE = ROOT / "All_auris_genes_beta_fungidb_rel70.tsv"
MASTER_FILE = ROOT / "candida_master_mapping.xlsx"
RESULTS_DIR = ROOT / "Results"
BACKUP_DIR = ROOT / "Backups"

MASTER_SHEET = "Master"

# FungiDB columns to add to the master workbook --> checked by Andy, use as final
COLUMNS_TO_ADD = {
    "Product Description": "FungiDB Product Description Rel70",
    "PFam ID": "PFam ID",
    "PFam Description": "PFam Description",
    "CDD ID": "CDD ID",
    "CDD Description": "CDD Description",
    "Gene3D ID": "Gene3D ID",
    "Gene3D Description": "Gene3D Description",
    "PANTHER ID": "PANTHER ID",
    "PANTHER Description": "PANTHER Description",
    "Prositefamilies ID": "ProSite ID",
    "Prositefamilies Description": "ProSite Description",
    "Smart ID": "SMART ID",
    "Smart Description": "SMART Description",
    "Superfamily ID": "Superfamily ID",
    "Superfamily Description": "Superfamily Description",
    "InterPro ID": "InterPro ID",
    "InterPro Description": "InterPro Description",
    "Computed GO Component IDs": "Computed GO Component IDs",
    "Computed GO Components": "Computed GO Components",
    "Computed GO Function IDs": "Computed GO Function IDs",
    "Computed GO Functions": "Computed GO Functions",
    "Computed GO Process IDs": "Computed GO Process IDs",
    "Computed GO Processes": "Computed GO Processes",
    "Curated GO Component IDs": "Curated GO Component IDs",
    "Curated GO Components": "Curated GO Components",
    "Curated GO Function IDs": "Curated GO Function IDs",
    "Curated GO Functions": "Curated GO Functions",
    "Curated GO Process IDs": "Curated GO Process IDs",
    "Curated GO Processes": "Curated GO Processes",
    "EC numbers": "EC numbers",
    "EC numbers from OrthoMCL": "EC numbers from OrthoMCL",
    "Ortholog count": "Ortholog count",
    "Paralog count": "Paralog count",
}

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


def clean_ids(series):
    return (
        series.fillna("")
        .astype(str)
        .str.strip()
        .replace("", pd.NA)
    )


def has_value(series):
    values = series.fillna("").astype(str).str.strip()

    return ~(
        values.eq("")
        | values.str.upper().eq("N/A")
        | values.str.casefold().eq("not assigned")
    )


def find_id_column(columns):
    candidates = [
        "Gene ID",
        "FungiDB Gene ID",
        "FungiDB ID",
        "VPDB Gene ID",
        "VPDB ID",
        "Beta FungiDB ID",
        "source_id",
    ]

    lower_columns = {
        str(column).strip().casefold(): column
        for column in columns
    }

    for candidate in candidates:
        if candidate.casefold() in lower_columns:
            return lower_columns[candidate.casefold()]

    return None


def main():
    for path in [FUNGIDB_FILE, MASTER_FILE]:
        if not path.exists():
            raise FileNotFoundError(path)

    RESULTS_DIR.mkdir(exist_ok=True)
    BACKUP_DIR.mkdir(exist_ok=True)

    # make a dated backup before changing the workbook
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    backup_file = BACKUP_DIR / f"candida_master_mapping_backup_{timestamp}.xlsx"
    shutil.copy2(MASTER_FILE, backup_file)

    fungidb = pd.read_csv(FUNGIDB_FILE, sep="\t", dtype=str)
    master = pd.read_excel(MASTER_FILE, sheet_name=MASTER_SHEET, dtype=str)

    if "Gene ID" not in fungidb.columns:
        raise KeyError("Gene ID column not found in FungiDB data")

    master_id = find_id_column(master.columns)
    if master_id is None:
        raise KeyError("Could not identify the gene ID column in the Master sheet")

    fungidb["Gene ID"] = clean_ids(fungidb["Gene ID"])
    master[master_id] = clean_ids(master[master_id])

    # classify hypothetical and annotation-poor proteins
    product = fungidb["Product Description"].fillna("").str.strip().str.casefold()
    fungidb["FungiDB Hypothetical"] = product.eq("hypothetical protein")

    found_evidence = [
        column for column in EVIDENCE_COLUMNS
        if column in fungidb.columns
    ]

    if not found_evidence:
        raise ValueError("No annotation evidence columns were found")

    evidence_matrix = pd.concat(
        [has_value(fungidb[column]) for column in found_evidence],
        axis=1,
    )
    fungidb["Has Conventional Evidence"] = evidence_matrix.any(axis=1)
    fungidb["Annotation Poor"] = (
        fungidb["FungiDB Hypothetical"]
        & ~fungidb["Has Conventional Evidence"]
    )

    # make sure the FungiDB IDs are unique before the merge
    fungidb_duplicates = fungidb["Gene ID"].dropna().duplicated().sum()
    master_duplicates = master[master_id].dropna().duplicated().sum()

    if fungidb_duplicates:
        raise ValueError("Duplicate FungiDB Gene IDs found")

    # keep the columns needed for the master workbook
    source_columns = [
        column for column in COLUMNS_TO_ADD
        if column in fungidb.columns
    ]

    lookup = fungidb[
        ["Gene ID", *source_columns]
    ].copy()

    lookup["FungiDB Hypothetical"] = fungidb["FungiDB Hypothetical"]
    lookup["Has Conventional Evidence"] = fungidb["Has Conventional Evidence"]
    lookup["Annotation Poor"] = fungidb["Annotation Poor"]

    lookup = lookup.rename(columns=COLUMNS_TO_ADD)
    lookup = lookup.rename(columns={"Gene ID": "__FungiDB_Match_ID"})

    # merge release 70 annotations into the master table
    merged = master.merge(
        lookup,
        how="left",
        left_on=master_id,
        right_on="__FungiDB_Match_ID",
        suffixes=("", "__NEW"),
        validate="m:1",
    )

    final_columns = list(COLUMNS_TO_ADD.values()) + [
        "FungiDB Hypothetical",
        "Has Conventional Evidence",
        "Annotation Poor",
    ]

    # prefer new release 70 values where the column already existed
    for column in final_columns:
        new_column = f"{column}__NEW"

        if new_column in merged.columns:
            merged[column] = merged[new_column].combine_first(merged[column])
            merged.drop(columns=new_column, inplace=True)

    merged.drop(columns="__FungiDB_Match_ID", inplace=True)

    # compare IDs and save anything that did not map
    fungidb_ids = set(fungidb["Gene ID"].dropna())
    master_ids = set(master[master_id].dropna())

    master_only = master_ids - fungidb_ids
    fungidb_only = fungidb_ids - master_ids

    if master_only:
        pd.DataFrame({
            "Master_ID_without_FungiDB_match": sorted(master_only)
        }).to_csv(
            RESULTS_DIR / "08_master_ids_without_fungidb_match.csv",
            index=False,
        )

    if fungidb_only:
        pd.DataFrame({
            "FungiDB_ID_without_Master_match": sorted(fungidb_only)
        }).to_csv(
            RESULTS_DIR / "08_fungidb_ids_without_master_match.csv",
            index=False,
        )

    # replace only the Master sheet and leave the other workbook sheets intact
    with pd.ExcelWriter(
        MASTER_FILE,
        engine="openpyxl",
        mode="a",
        if_sheet_exists="replace",
    ) as writer:
        merged.to_excel(writer, sheet_name=MASTER_SHEET, index=False)

    print(f"hypothetical proteins: {int(fungidb['FungiDB Hypothetical'].sum())}")
    print(
        "hypothetical proteins with evidence: "
        f"{int((fungidb['FungiDB Hypothetical'] & fungidb['Has Conventional Evidence']).sum())}"
    )
    print(f"annotation-poor proteins: {int(fungidb['Annotation Poor'].sum())}")
    print(f"duplicate FungiDB IDs: {fungidb_duplicates}")
    print(f"duplicate Master IDs: {master_duplicates}")
    print(f"matched IDs: {len(master_ids & fungidb_ids)}")
    print(f"Master-only IDs: {len(master_only)}")
    print(f"FungiDB-only IDs: {len(fungidb_only)}")


if __name__ == "__main__":
    main()


# GAI declaration
# GAI was used in this script for lines 98-118.
# i asked for a generated response on how to find an identifier column when the
# same type of ID may have different column names in different files. i was given
# a generic solution using a list such as "your_possible_column_names", making
# a case-insensitive dictionary from "your_dataframe.columns", then returning
# the first matching name. that was taken and adapted into find_id_column() to
# recognise the different B8441 ID column names used in this project.

# GAI was used in this script for lines 151-167.
# i asked for a generated response on how to check a list of annotation columns
# and create one boolean column showing whether any evidence is present for each
# row. i was given a generic solution which applied a function to each column,
# combined the boolean Series using pd.concat(..., axis=1), then used
# .any(axis=1). that was taken and adapted to create Has Conventional Evidence
# from the available FungiDB annotation sources.

# GAI was used in this script for lines 193-217.
# i asked for a generated response on how to merge new annotation data into an
# existing dataframe when some column names already exist, while keeping the new
# value where available and falling back to the old value where it is missing.
# i was given a generic solution using merge(..., suffixes=("", "__NEW")) and
# "new_column".combine_first("old_column"). that was taken and adapted to update
# the Release 70 annotation columns in the Master mapping.

# GAI was used in this script for lines 219-240.
# i asked for a generated response on how to find identifiers that occur in one
# dataframe but not another. i was given a generic solution converting
# "dataframe_1['your_id']" and "dataframe_2['your_id']" into sets and using set
# subtraction to find unmatched values. that was taken and adapted to create the
# Master-only and FungiDB-only ID result files.

# GAI was used in this script for lines 242-249.
# i asked for a generated response on how to replace one worksheet in an existing
# Excel workbook from a pandas dataframe without overwriting the other sheets.
# i was given a generic example using pd.ExcelWriter("your_file_path",
# engine="openpyxl", mode="a", if_sheet_exists="replace") followed by
# "your_dataframe".to_excel(). that was taken and adapted to replace only the
# Master sheet in candida_master_mapping.xlsx.