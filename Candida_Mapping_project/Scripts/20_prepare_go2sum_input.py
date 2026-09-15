from pathlib import Path

import pandas as pd

# file paths
ROOT = Path(__file__).resolve().parents[1]

PROTNOTE_FILE = ROOT / "Results" / "12_protnote_annotation_poor_all_predictions.csv"
GO2SUM_FILE = ROOT / "Results" / "20_go2sum_protnote_annotation_poor_input.tab"
TEST_FILE = ROOT / "Results" / "20_go2sum_test_5_proteins.tab"
UNIPROT_TEST_FILE = ROOT / "Results" / "20_go2sum_test_uniprot_id.tab"


def combine_go_ids(values):
    # keep GO IDs in their original order without duplicates
    cleaned = [
        str(value).strip()
        for value in values
        if pd.notna(value) and str(value).strip()
    ]
    return ";".join(dict.fromkeys(cleaned))


def main():
    if not PROTNOTE_FILE.exists():
        raise FileNotFoundError(PROTNOTE_FILE)

    df = pd.read_csv(PROTNOTE_FILE)

    needed = {"sequence_id", "go_id", "Is_Obsolete"}
    missing = needed - set(df.columns)

    if missing:
        raise KeyError(f"Missing ProtNote columns: {', '.join(sorted(missing))}")

    # keep GO predictions that are present and not marked obsolete
    obsolete = (
        df["Is_Obsolete"]
        .fillna(False)
        .astype(str)
        .str.casefold()
        .eq("true")
    )

    usable = df[
        df["go_id"].notna()
        & df["go_id"].astype(str).str.strip().ne("")
        & ~obsolete
    ].copy()

    # GO2SUM expects one protein per row with semicolon-separated GO IDs
    go2sum = (
        usable
        .groupby("sequence_id")["go_id"]
        .apply(combine_go_ids)
        .reset_index(name="GO_IDs")
        .rename(columns={"sequence_id": "Protein"})
    )

    GO2SUM_FILE.parent.mkdir(exist_ok=True)
    go2sum.to_csv(GO2SUM_FILE, sep="\t", index=False)

    # small batch used while troubleshooting the GO2SUM server
    go2sum.head(5).to_csv(TEST_FILE, sep="\t", index=False)

    # diagnostic file using a recognised UniProt-style identifier
    uniprot_test = go2sum.head(1).copy()
    if not uniprot_test.empty:
        uniprot_test.loc[:, "Protein"] = "Q5AK66"
    uniprot_test.to_csv(UNIPROT_TEST_FILE, sep="\t", index=False)

    term_counts = go2sum["GO_IDs"].str.count(";").add(1)

    print(f"ProtNote rows: {len(df)}")
    print(f"unique proteins: {df['sequence_id'].nunique()}")
    print(f"non-obsolete GO predictions: {len(usable)}")
    print(f"GO2SUM proteins: {len(go2sum)}")

    if not term_counts.empty:
        print(f"mean GO terms per protein: {term_counts.mean():.2f}")
        print(f"median GO terms per protein: {term_counts.median():.1f}")


if __name__ == "__main__":
    main()


# GAI declaration
# GAI was used for lines 14-21 when combining repeated GO predictions into one
# semicolon-separated field without losing their original order. i asked how to
# remove duplicates while preserving order, and the example used:
#
# unique = dict.fromkeys(your_values)
# combined = ";".join(unique)
#
# i adapted this into combine_go_ids() after first removing blank values.

# i also used GAI around lines 52-58 to reshape the prediction table into one row
# per protein. the suggested pandas pattern was:
#
# result = df.groupby("your_id")["your_value"].apply(your_function).reset_index()
#
# this was adapted to sequence_id and go_id to create the GO2SUM input table.