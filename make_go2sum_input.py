#!/usr/bin/env python3

"""
Convert predicted GO terms into a simple GO2Sum-style input table.

Input TSV columns expected:
protein_id, protnote_go

Example:
python make_go2sum_input.py example_data/tiny_cauris_example.tsv go2sum_input.tsv
"""

import sys
import pandas as pd


def main():
    if len(sys.argv) != 3:
        print("Usage: python make_go2sum_input.py input.tsv output.tsv")
        sys.exit(1)

    input_file = sys.argv[1]
    output_file = sys.argv[2]

    df = pd.read_csv(input_file, sep="\t")

    if "protein_id" not in df.columns or "protnote_go" not in df.columns:
        raise ValueError("Input file must contain protein_id and protnote_go columns")

    rows = []

    for _, row in df.iterrows():
        protein_id = row["protein_id"]
        go_terms = str(row["protnote_go"]).strip()

        if go_terms == "" or go_terms.lower() == "nan":
            continue

        rows.append({
            "protein_id": protein_id,
            "go_terms": go_terms.replace(";", ",")
        })

    out = pd.DataFrame(rows)
    out.to_csv(output_file, sep="\t", index=False)

    print(f"Saved GO2Sum input for {len(out)} proteins to {output_file}")


if __name__ == "__main__":
    main()