from pathlib import Path

import pandas as pd

# file paths
ROOT = Path(__file__).resolve().parents[1]

PROTNOTE_FILE = ROOT / "Data" / "raw" / "protnote" / "protnote_annotation_poor_predictions.tsv"
CANDIDATE_FILE = ROOT / "Results" / "fungidb_rel70_annotation_poor_candidate_list.csv"
RESULTS_DIR = ROOT / "Results"

OUTPUTS = {
    "all": RESULTS_DIR / "12_protnote_annotation_poor_all_predictions.csv",
    "specific": RESULTS_DIR / "12_protnote_specific_predictions.csv",
    "top": RESULTS_DIR / "12_protnote_top_specific_prediction_per_protein.csv",
    "high_conf": RESULTS_DIR / "12_protnote_high_confidence_specific_predictions.csv",
    "missing": RESULTS_DIR / "12_annotation_poor_candidates_without_protnote.csv",
    "candidates": RESULTS_DIR / "12_protnote_candidate_predictions_for_validation.csv",
}

GO_COLUMNS = [
    "Computed GO Component IDs",
    "Computed GO Function IDs",
    "Computed GO Process IDs",
]

GENERIC_GO_IDS = {
    "GO:0008150",
    "GO:0003674",
    "GO:0005575",
    "GO:0044464",
    "GO:0044425",
    "GO:0043226",
    "GO:0043227",
    "GO:0044444",
    "GO:0005488",
    "GO:0003824",
    "GO:0009987",
    "GO:0065007",
}

BROAD_DESCRIPTIONS = {
    "intracellular part",
    "intracellular organelle part",
    "organelle part",
    "intracellular organelle",
    "intracellular membrane-bounded organelle",
    "membrane",
    "metabolic process",
    "organic substance metabolic process",
    "cellular metabolic process",
    "nitrogen compound metabolic process",
    "macromolecule metabolic process",
    "primary metabolic process",
    "regulation of biological process",
    "regulation of cellular process",
    "cellular component organization",
    "cellular component organization or biogenesis",
    "protein-containing complex",
    "localization",
    "cytoplasm",
    "nucleus",
}

HIGH_CONF = 0.8


def has_annotation(series):
    values = series.fillna("").astype(str).str.strip()

    return ~(
        values.eq("")
        | values.str.upper().eq("N/A")
        | values.str.casefold().isin(["not assigned", "nan"])
    )


def add_context(data, context):
    return data.merge(
        context,
        how="left",
        left_on="sequence_id",
        right_on="Gene ID",
    )


def main():
    for path in [PROTNOTE_FILE, CANDIDATE_FILE]:
        if not path.exists():
            raise FileNotFoundError(path)

    RESULTS_DIR.mkdir(exist_ok=True)

    protnote = pd.read_csv(PROTNOTE_FILE, sep="\t")
    annotation = pd.read_csv(CANDIDATE_FILE, dtype=str)

    needed = {"sequence_id", "go_id", "probability", "description", "rank"}
    missing = needed - set(protnote.columns)

    if missing:
        raise KeyError(f"Missing ProtNote columns: {', '.join(sorted(missing))}")

    if "Gene ID" not in annotation.columns:
        raise KeyError("Gene ID column not found in candidate file")

    protnote["probability"] = pd.to_numeric(
        protnote["probability"],
        errors="coerce",
    )
    protnote["rank"] = pd.to_numeric(
        protnote["rank"],
        errors="coerce",
    )

    # check ProtNote IDs against the annotation-poor candidate set
    protnote_ids = set(
        protnote["sequence_id"]
        .dropna()
        .astype(str)
        .str.strip()
    )
    candidate_ids = set(
        annotation["Gene ID"]
        .dropna()
        .astype(str)
        .str.strip()
    )

    matched = protnote_ids & candidate_ids
    not_predicted = candidate_ids - protnote_ids
    unexpected = protnote_ids - candidate_ids

    pd.DataFrame({"Gene ID": sorted(not_predicted)}).to_csv(
        OUTPUTS["missing"],
        index=False,
    )

    # confirm that the candidate set still lacks computed GO annotation
    existing_go = pd.Series(False, index=annotation.index)

    for column in GO_COLUMNS:
        if column in annotation.columns:
            existing_go |= has_annotation(annotation[column])

    # remove very broad/root terms and obsolete GO predictions
    protnote["Is_Generic"] = (
        protnote["go_id"]
        .astype(str)
        .isin(GENERIC_GO_IDS)
    )

    protnote["Is_Obsolete"] = (
        protnote["description"]
        .fillna("")
        .astype(str)
        .str.casefold()
        .str.contains("obsolete")
    )

    protnote["Is_Specific"] = (
        ~protnote["Is_Generic"]
        & ~protnote["Is_Obsolete"]
    )

    specific = protnote[protnote["Is_Specific"]].copy()

    # keep the highest-ranked specific prediction for each protein
    top_specific = (
        specific
        .sort_values(
            ["sequence_id", "rank", "probability"],
            ascending=[True, True, False],
        )
        .groupby("sequence_id", as_index=False)
        .first()
    )

    high_conf = specific[
        specific["probability"] >= HIGH_CONF
    ].copy()

    # add useful protein information from the candidate table
    context_columns = [
        "Gene ID",
        "source_id",
        "Product Description",
        "Protein Length",
        "Ortholog count",
        "Paralog count",
    ]
    context_columns = [
        column
        for column in context_columns
        if column in annotation.columns
    ]
    context = annotation[context_columns].copy()

    all_with_context = add_context(protnote, context)
    specific_with_context = add_context(specific, context)
    top_with_context = add_context(top_specific, context)
    high_conf_with_context = add_context(high_conf, context)

    # remove broad terms from the high-confidence validation candidates
    descriptions = (
        high_conf_with_context["description"]
        .fillna("")
        .str.strip()
        .str.casefold()
    )

    validation_candidates = high_conf_with_context[
        ~descriptions.isin(BROAD_DESCRIPTIONS)
    ].copy()

    validation_candidates = validation_candidates.sort_values(
        ["probability", "rank"],
        ascending=[False, True],
    )

    all_with_context.to_csv(OUTPUTS["all"], index=False)
    specific_with_context.to_csv(OUTPUTS["specific"], index=False)
    top_with_context.to_csv(OUTPUTS["top"], index=False)
    high_conf_with_context.to_csv(OUTPUTS["high_conf"], index=False)
    validation_candidates.to_csv(OUTPUTS["candidates"], index=False)

    predictions_per_protein = protnote.groupby("sequence_id").size()

    print(f"ProtNote rows: {len(protnote)}")
    print(f"proteins with predictions: {protnote['sequence_id'].nunique()}")
    print(f"matched candidate proteins: {len(matched)}")
    print(f"candidates without ProtNote output: {len(not_predicted)}")
    print(f"unexpected ProtNote proteins: {len(unexpected)}")
    print(f"candidates with existing computed GO: {int(existing_go.sum())}")
    print(
        f"predictions per protein: "
        f"{predictions_per_protein.min()}-{predictions_per_protein.max()}"
    )
    print(f"generic predictions: {int(protnote['Is_Generic'].sum())}")
    print(f"obsolete predictions: {int(protnote['Is_Obsolete'].sum())}")
    print(f"specific predictions: {len(specific)}")
    print(
        f"proteins with a specific prediction: "
        f"{specific['sequence_id'].nunique()}"
    )
    print(f"high-confidence specific predictions: {len(high_conf)}")
    print(f"validation candidates: {len(validation_candidates)}")


if __name__ == "__main__":
    main()


# GAI declaration
# for lines 121-136, i asked ChatGPT how to compare two groups of protein IDs
# and separate them into matched IDs, expected proteins with no prediction, and
# unexpected predictions. the example used Python sets called "expected_ids"
# and "predicted_ids" with intersection and subtraction. i adapted that to the
# candidate Gene IDs and ProtNote sequence IDs used here.

# GAI was used for lines 168-178 when deciding how to retain one specific
# ProtNote result per protein. i asked how to sort multiple predictions by ID,
# rank and confidence, then select the first row for each protein. the response
# showed a general pandas example using sort_values() followed by
# groupby("your_id", as_index=False).first(), which was adapted to sequence_id,
# rank and probability.

# i asked another GAI question for lines 207-223 about removing broad GO
# descriptions from the high-confidence results. rather than writing a long
# series of comparisons, the suggested example cleaned "your_text_column" with
# strip/lower and filtered using ~.isin(your_excluded_terms). this was applied
# to the BROAD_DESCRIPTIONS set before producing the validation candidate table.

# the helper on lines 78-84 and its use on lines 198-201 also came from GAI
# assistance. i asked how to avoid repeating the same merge several times when
# adding the same protein metadata to different result tables. ChatGPT suggested
# wrapping the merge in a small function taking "your_data" and "your_context",
# which became add_context() and is reused for each ProtNote subset.