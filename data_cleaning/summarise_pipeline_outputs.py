#!/usr/bin/env python3

"""
Summarise exploratory annotation results.

Example:
python summarise_pipeline_outputs.py comparison_summary.tsv
"""

import sys
import pandas as pd


def main():
    if len(sys.argv) != 2:
        print("Usage: python summarise_pipeline_outputs.py comparison_summary.tsv")
        sys.exit(1)

    df = pd.read_csv(sys.argv[1], sep="\t")

    print("=== C. auris annotation prototype summary ===")
    print(f"Total proteins checked: {len(df)}")

    print("\nStatus counts:")
    print(df["status"].value_counts())

    rescue = df[df["status"] == "high_priority_rescue_candidate"]
    print(f"\nHigh-priority rescue candidates: {len(rescue)}")

    if len(rescue) > 0:
        print("\nExample rescue candidates:")
        print(rescue[["protein_id", "current_name", "predicted_go_count", "has_cath_or_foldseek"]].head(10).to_string(index=False))


if __name__ == "__main__":
    main()