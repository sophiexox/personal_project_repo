from pathlib import Path

import pandas as pd

# file paths
ROOT = Path(__file__).resolve().parents[1]

PROTNOTE_FILE = ROOT / "Results" / "12_protnote_annotation_poor_all_predictions.csv"
PROTNLM_FILE = ROOT / "Data" / "raw" / "protnlm" / "ProtNLM_annotation_poor_521_summary.csv"
OUTPUT_FILE = ROOT / "Results" / "21_integrated_annotation_poor_predictions.csv"

HIGH_CONF_THRESHOLD = 0.5

GENERIC_TERMS = [
    "hypothetical protein",
    "hypothetical_protein",
    "conserved protein",
    "conserved hypothetical protein",
    "expressed protein",
    "similar to",
    "unknown protein",
    "uncharacterized protein",
    "uncharacterised protein",
]


def as_bool(series):
    # handle boolean columns after they have been written to and read from CSV
    return (
        series.fillna(False)
        .astype(str)
        .str.strip()
        .str.casefold()
        .eq("true")
    )


def is_generic(value):
    if pd.isna(value):
        return True

    text = str(value).strip().casefold()
    return any(
        text == term.casefold() or text.startswith(term.casefold())
        for term in GENERIC_TERMS
    )


def top_prediction(data, names):
    # keep the highest-ranked prediction for each protein
    top = (
        data
        .sort_values(["sequence_id", "rank"])
        .groupby("sequence_id", as_index=False)
        .first()
    )

    return top[
        ["sequence_id", "go_id", "description", "probability", "rank"]
    ].rename(columns={
        "sequence_id": "Protein",
        "go_id": names["go"],
        "description": names["description"],
        "probability": names["probability"],
        "rank": names["rank"],
    })


def main():
    for path in [PROTNOTE_FILE, PROTNLM_FILE]:
        if not path.exists():
            raise FileNotFoundError(path)

    protnote = pd.read_csv(PROTNOTE_FILE)
    protnlm = pd.read_csv(PROTNLM_FILE)

    needed_protnote = {
        "sequence_id",
        "go_id",
        "description",
        "probability",
        "rank",
        "Is_Obsolete",
        "Is_Specific",
    }
    missing = needed_protnote - set(protnote.columns)

    if missing:
        raise KeyError(f"Missing ProtNote columns: {', '.join(sorted(missing))}")

    needed_protnlm = {"Protein", "ProtNLM_Top1", "ProtNLM_Top1_Score"}
    missing = needed_protnlm - set(protnlm.columns)

    if missing:
        raise KeyError(f"Missing ProtNLM columns: {', '.join(sorted(missing))}")

    protnote["probability"] = pd.to_numeric(
        protnote["probability"],
        errors="coerce",
    )
    protnote["rank"] = pd.to_numeric(protnote["rank"], errors="coerce")
    protnlm["ProtNLM_Top1_Score"] = pd.to_numeric(
        protnlm["ProtNLM_Top1_Score"],
        errors="coerce",
    )

    protnote["Is_Obsolete"] = as_bool(protnote["Is_Obsolete"])
    protnote["Is_Specific"] = as_bool(protnote["Is_Specific"])

    # summarise usable ProtNote predictions to one row per protein
    usable = protnote[
        protnote["go_id"].notna()
        & ~protnote["Is_Obsolete"]
    ].copy()
    specific = usable[usable["Is_Specific"]].copy()

    protnote_top = top_prediction(
        usable,
        {
            "go": "ProtNote_Top_GO",
            "description": "ProtNote_Top_GO_Description",
            "probability": "ProtNote_Top_Probability",
            "rank": "ProtNote_Top_Rank",
        },
    )
    protnote_top_specific = top_prediction(
        specific,
        {
            "go": "ProtNote_Top_Specific_GO",
            "description": "ProtNote_Top_Specific_Description",
            "probability": "ProtNote_Top_Specific_Probability",
            "rank": "ProtNote_Top_Specific_Rank",
        },
    )

    counts = (
        usable
        .groupby("sequence_id")
        .agg(
            ProtNote_NonObsolete_Count=("go_id", "count"),
            ProtNote_Specific_Count=("Is_Specific", "sum"),
        )
        .reset_index()
        .rename(columns={"sequence_id": "Protein"})
    )

    high_conf_counts = (
        specific[specific["probability"] >= HIGH_CONF_THRESHOLD]
        .groupby("sequence_id")
        .size()
        .reset_index(name="ProtNote_HighConf_Specific_Count")
        .rename(columns={"sequence_id": "Protein"})
    )

    specific_descriptions = (
        specific
        .sort_values(["sequence_id", "rank"])
        .groupby("sequence_id")["description"]
        .apply(
            lambda values: " | ".join(
                dict.fromkeys(values.dropna().astype(str))
            )
        )
        .reset_index(name="ProtNote_All_Specific_Descriptions")
        .rename(columns={"sequence_id": "Protein"})
    )

    # recover protein metadata already carried through the ProtNote result file
    metadata_columns = [
        column
        for column in [
            "sequence_id",
            "Gene ID",
            "source_id",
            "Product Description",
            "Protein Length",
            "Ortholog count",
            "Paralog count",
        ]
        if column in protnote.columns
    ]

    metadata = (
        protnote[metadata_columns]
        .drop_duplicates(subset="sequence_id")
        .rename(columns={"sequence_id": "Protein"})
    )

    # add ProtNLM genericity and score flags used only for candidate ranking
    protnlm["ProtNLM_Top1_Generic"] = protnlm["ProtNLM_Top1"].apply(is_generic)
    protnlm["ProtNLM_Score_GE_0.2"] = protnlm["ProtNLM_Top1_Score"] >= 0.2
    protnlm["ProtNLM_Score_GE_0.5"] = protnlm["ProtNLM_Top1_Score"] >= 0.5
    protnlm["ProtNLM_Score_GE_0.8"] = protnlm["ProtNLM_Top1_Score"] >= 0.8

    integrated = metadata.merge(protnlm, on="Protein", how="outer")

    for table in [
        protnote_top,
        protnote_top_specific,
        counts,
        high_conf_counts,
        specific_descriptions,
    ]:
        integrated = integrated.merge(table, on="Protein", how="left")

    count_columns = [
        "ProtNote_NonObsolete_Count",
        "ProtNote_Specific_Count",
        "ProtNote_HighConf_Specific_Count",
    ]
    integrated[count_columns] = (
        integrated[count_columns]
        .fillna(0)
        .astype(int)
    )

    # ranking flags for manual candidate inspection, not biological confidence
    integrated["Has_ProtNote_Specific_Prediction"] = (
        integrated["ProtNote_Specific_Count"] > 0
    )
    integrated["Has_HighConf_ProtNote_Specific"] = (
        integrated["ProtNote_HighConf_Specific_Count"] > 0
    )
    integrated["ProtNLM_NonGeneric"] = ~(
        integrated["ProtNLM_Top1_Generic"]
        .fillna(True)
        .astype(bool)
    )

    score_flags = [
        "ProtNLM_Score_GE_0.2",
        "ProtNLM_Score_GE_0.5",
        "ProtNLM_NonGeneric",
        "Has_ProtNote_Specific_Prediction",
        "Has_HighConf_ProtNote_Specific",
    ]
    integrated["Candidate_Priority_Score"] = (
        integrated[score_flags]
        .fillna(False)
        .astype(int)
        .sum(axis=1)
    )

    integrated = integrated.sort_values(
        [
            "Candidate_Priority_Score",
            "ProtNLM_Top1_Score",
            "ProtNote_Top_Specific_Probability",
        ],
        ascending=[False, False, False],
    )

    OUTPUT_FILE.parent.mkdir(exist_ok=True)
    integrated.to_csv(OUTPUT_FILE, index=False)

    print(f"integrated proteins: {integrated['Protein'].nunique()}")
    print(f"proteins with ProtNLM: {integrated['ProtNLM_Top1'].notna().sum()}")
    print(
        "proteins with specific ProtNote prediction: "
        f"{integrated['Has_ProtNote_Specific_Prediction'].sum()}"
    )
    print(
        "proteins with high-confidence specific ProtNote prediction: "
        f"{integrated['Has_HighConf_ProtNote_Specific'].sum()}"
    )
    print(f"non-generic ProtNLM top-1: {integrated['ProtNLM_NonGeneric'].sum()}")

    print("\npriority score distribution:")
    print(
        integrated["Candidate_Priority_Score"]
        .value_counts()
        .sort_index()
        .to_string()
    )


if __name__ == "__main__":
    main()


# GAI declaration
# GAI was used for lines 49-66 to avoid repeating the same sorting, grouping and
# renaming code for the two ProtNote summaries. i asked how one helper could take
# a dataframe and a set of output column names, then return the first ranked row
# per protein. the example used sort_values(), groupby().first() and rename(),
# which was adapted into top_prediction().

# lines 230-242 were also developed with GAI support. i wanted the candidate
# priority score to be the number of True ranking flags rather than five separate
# `.loc[...] += 1` statements. the example response converted several boolean
# columns to integers and summed across rows:
#
# df["score"] = (
#     df[your_flag_columns]
#     .fillna(False)
#     .astype(int)
#     .sum(axis=1)
# )
#
# i adapted that to the ProtNLM and ProtNote flags used for manual prioritisation.