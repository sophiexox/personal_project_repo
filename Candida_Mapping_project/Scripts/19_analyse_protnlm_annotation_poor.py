from pathlib import Path

import pandas as pd

# file paths
ROOT = Path(__file__).resolve().parents[1]

PROTNLM_DIR = ROOT / "Data" / "raw" / "protnlm"
PREDICTIONS_FILE = PROTNLM_DIR / "ProtNLM_annotation_poor_521_predictions.csv"
SUMMARY_FILE = PROTNLM_DIR / "ProtNLM_annotation_poor_521_summary.csv"
OUTPUT_FILE = ROOT / "Results" / "19_protnlm_annotation_poor_summary_analysis.csv"

SCORE_THRESHOLDS = [0.9, 0.8, 0.5, 0.2, 0.1]

GENERIC_TERMS = [
    "hypothetical protein",
    "hypothetical_protein",
    "uncharacterized protein",
    "uncharacterised protein",
    "expressed protein",
    "conserved protein",
    "similar to",
]


def is_generic(label):
    # flag broad labels that do not give useful functional information
    text = str(label).strip().casefold()
    return any(term in text for term in GENERIC_TERMS)


def main():
    for path in [PREDICTIONS_FILE, SUMMARY_FILE]:
        if not path.exists():
            raise FileNotFoundError(path)

    predictions = pd.read_csv(PREDICTIONS_FILE)
    summary = pd.read_csv(SUMMARY_FILE)

    required = {"Protein", "ProtNLM_Top1", "ProtNLM_Top1_Score"}
    missing = required - set(summary.columns)

    if missing:
        raise KeyError(f"Missing summary columns: {', '.join(sorted(missing))}")

    summary["ProtNLM_Top1_Score"] = pd.to_numeric(
        summary["ProtNLM_Top1_Score"],
        errors="raise",
    )

    # normalise labels and flag generic top predictions
    summary["ProtNLM_Top1_normalised"] = (
        summary["ProtNLM_Top1"]
        .fillna("")
        .astype(str)
        .str.strip()
        .str.casefold()
    )
    summary["ProtNLM_Generic"] = summary["ProtNLM_Top1"].apply(is_generic)

    OUTPUT_FILE.parent.mkdir(exist_ok=True)
    summary.to_csv(OUTPUT_FILE, index=False)

    print(f"prediction rows: {len(predictions)}")
    print(f"proteins in full predictions: {predictions['Protein'].nunique()}")
    print(f"proteins in summary: {summary['Protein'].nunique()}")
    print(f"mean top-1 score: {summary['ProtNLM_Top1_Score'].mean():.4f}")
    print(f"median top-1 score: {summary['ProtNLM_Top1_Score'].median():.4f}")

    print("\nscore thresholds:")
    for threshold in SCORE_THRESHOLDS:
        count = int((summary["ProtNLM_Top1_Score"] >= threshold).sum())
        percentage = count / len(summary) * 100
        print(f">= {threshold:.1f}: {count} ({percentage:.1f}%)")

    generic_count = int(summary["ProtNLM_Generic"].sum())
    generic_pct = generic_count / len(summary) * 100
    print(f"\ngeneric top-1 predictions: {generic_count} ({generic_pct:.1f}%)")

    print("\nmost common top-1 predictions:")
    print(summary["ProtNLM_Top1"].value_counts().head(25).to_string())

    print("\nhighest-scoring top-1 predictions:")
    print(
        summary[
            ["Protein", "ProtNLM_Top1", "ProtNLM_Top1_Score", "ProtNLM_Generic"]
        ]
        .sort_values("ProtNLM_Top1_Score", ascending=False)
        .head(20)
        .to_string(index=False)
    )


if __name__ == "__main__":
    main()


# GAI declaration
# GAI was used for lines 26-29 when making the generic-label check reusable.
# i asked how to test whether any phrase from a list occurred inside a protein
# label without writing a separate comparison for every phrase. the example used:
#
# text = str(your_label).lower()
# return any(term in text for term in your_terms)
#
# i adapted this into is_generic() for the ProtNLM top-1 labels.

# lines 70-74 also use GAI-supported code. i asked how to run the same count
# across several score cut-offs and calculate the corresponding percentage:
#
# for threshold in thresholds:
#     count = (df["score"] >= threshold).sum()
#     percent = count / len(df) * 100
#
# this was adapted to ProtNLM_Top1_Score and the thresholds used in the analysis.