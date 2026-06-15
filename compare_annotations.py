#!/usr/bin/env python3

"""
Compare existing C. auris annotations against predicted GO-based annotations.

This is an early validation scaffold:
- Finds proteins where InterPro and ProtNote/GO predictions agree.
- Finds proteins rescued by predicted GO terms when current annotation is weak.
- Flags cases where CATH/FoldSeek and GO prediction support the same broad function.

Example:
python compare_annotations.py example_data/tiny_cauris_example.tsv comparison_summary.tsv
"""

import sys
import pandas as pd


def split_terms(value):
    if pd.isna(value) or str(value).strip() == "":
        return set()
    return set(term.strip() for term in str(value).replace(",", ";").split(";") if term.strip())


def is_weak_name(name):
    weak_words = ["hypothetical", "uncharacterised", "uncharacterized", "unknown", "predicted protein"]
    name = str(name).lower()
    return any(word in name for word in weak_words)


def main():
    if len(sys.argv) != 3:
        print("Usage: python compare_annotations.py input.tsv output.tsv")
        sys.exit(1)

    input_file = sys.argv[1]
    output_file = sys.argv[2]

    df = pd.read_csv(input_file, sep="\t")

    results = []

    for _, row in df.iterrows():
        protein_id = row["protein_id"]
        current_name = row.get("current_name", "")

        interpro_go = split_terms(row.get("interpro_go", ""))
        predicted_go = split_terms(row.get("protnote_go", ""))

        overlap = interpro_go.intersection(predicted_go)

        weak_current_annotation = is_weak_name(current_name)
        has_predicted_go = len(predicted_go) > 0
        has_cath = str(row.get("cath_label", "")).strip() not in ["", "nan"]

        if weak_current_annotation and has_predicted_go and has_cath:
            status = "high_priority_rescue_candidate"
        elif len(overlap) > 0:
            status = "supported_by_existing_GO"
        elif has_predicted_go:
            status = "new_prediction_needs_validation"
        else:
            status = "no_prediction"

        results.append({
            "protein_id": protein_id,
            "current_name": current_name,
            "interpro_go_count": len(interpro_go),
            "predicted_go_count": len(predicted_go),
            "go_overlap_count": len(overlap),
            "has_cath_or_foldseek": has_cath,
            "status": status,
            "overlapping_go_terms": ";".join(sorted(overlap))
        })

    out = pd.DataFrame(results)
    out.to_csv(output_file, sep="\t", index=False)

    print(out["status"].value_counts())
    print(f"\nSaved comparison table to {output_file}")


if __name__ == "__main__":
    main()