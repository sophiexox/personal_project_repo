from pathlib import Path

import pandas as pd

# file paths
ROOT = Path(__file__).resolve().parents[1]

INPUT_FILE = ROOT / "Results" / "12_protnote_candidate_predictions_for_validation.csv"
OUTPUT_FILE = ROOT / "Results" / "13_protnote_candidate_summary_by_protein.csv"

SUSPICIOUS_WORDS = [
    "viral",
    "virus",
    "immune",
    "host",
    "pathogenesis",
    "virulence",
    "development",
    "cell population proliferation",
]

FUNCTION_WORDS = [
    "dna",
    "rna",
    "nucleic acid",
    "transcription",
    "translation",
    "kinase",
    "phosphatase",
    "atpase",
    "transport",
    "transporter",
    "oxidoreductase",
    "hydrolase",
    "protease",
    "peptidase",
    "metal ion",
    "mitochond",
    "cell wall",
]


def contains_keyword(text, keywords):
    text = str(text).casefold()
    return any(word.casefold() in text for word in keywords)


def main():
    if not INPUT_FILE.exists():
        raise FileNotFoundError(INPUT_FILE)

    df = pd.read_csv(INPUT_FILE)

    needed = {"sequence_id", "probability", "rank", "description", "go_id"}
    missing = needed - set(df.columns)

    if missing:
        raise KeyError(f"Missing columns: {', '.join(sorted(missing))}")

    df["probability"] = pd.to_numeric(df["probability"], errors="coerce")
    df["rank"] = pd.to_numeric(df["rank"], errors="coerce")

    # make one summary row for each protein
    summary_rows = []

    for sequence_id, group in df.groupby("sequence_id"):
        group = group.sort_values(
            ["probability", "rank"],
            ascending=[False, True],
        )

        top = group.head(5)
        first = top.iloc[0]

        prediction_summary = " | ".join(
            f"{row.description} ({row.probability:.4f})"
            for row in top.itertuples()
        )
        go_summary = " | ".join(top["go_id"].astype(str))

        summary_rows.append({
            "Gene ID": sequence_id,
            "Protein Length": first.get("Protein Length"),
            "Highest Probability": group["probability"].max(),
            "Number High-Confidence Predictions": len(group),
            "Top Candidate GO ID": first["go_id"],
            "Top Candidate Description": first["description"],
            "Top Candidate Probability": first["probability"],
            "Top 5 Candidate Predictions": prediction_summary,
            "Top 5 GO IDs": go_summary,
        })

    summary = pd.DataFrame(summary_rows)

    # simple keyword flags for manual candidate review
    summary["Flag_Unusual_or_Context_Sensitive"] = (
        summary["Top 5 Candidate Predictions"]
        .apply(lambda text: contains_keyword(text, SUSPICIOUS_WORDS))
    )
    summary["Flag_Specific_Function"] = (
        summary["Top 5 Candidate Predictions"]
        .apply(lambda text: contains_keyword(text, FUNCTION_WORDS))
    )

    summary = summary.sort_values(
        [
            "Flag_Unusual_or_Context_Sensitive",
            "Flag_Specific_Function",
            "Highest Probability",
        ],
        ascending=[False, False, False],
    )

    summary.to_csv(OUTPUT_FILE, index=False)

    print(f"proteins summarised: {len(summary)}")
    print(
        "unusual/context-sensitive flags: "
        f"{int(summary['Flag_Unusual_or_Context_Sensitive'].sum())}"
    )
    print(
        "specific-function flags: "
        f"{int(summary['Flag_Specific_Function'].sum())}"
    )


if __name__ == "__main__":
    main()


# GAI declaration
# GAI was used for lines 63-91. i asked how to reduce multiple prediction rows
# for each protein to a single summary while keeping the five strongest results.
# ChatGPT showed a general groupby("your_id") loop, sorting each group by a score
# column and using head(5). i adapted this to sequence_id, probability and rank,
# then kept the top GO term plus a combined top-five summary.

# for lines 75-79, i asked how to combine text and numeric values from several
# dataframe rows into one readable string. the example response used
# " | ".join(f"{row.your_text} ({row.your_score:.4f})" for row in
# your_dataframe.itertuples()). i used the same idea for the ProtNote
# descriptions and probabilities.

# lines 43-45 and 95-103 use another GAI suggestion. i needed a quick way to flag
# a protein when any word from a keyword list appeared in its combined prediction
# text. the suggested generic solution used any(keyword in text for keyword in
# keywords), with both sides converted to lower case. i adapted this into
# contains_keyword() and applied it to the two candidate-review keyword lists.