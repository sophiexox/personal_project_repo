from pathlib import Path

import pandas as pd

# file paths
ROOT = Path(__file__).resolve().parents[1]

INPUT_FILE = ROOT / "Data" / "raw" / "uniprot" / "uniprot_candidozyma_auris.tsv"
OUTPUT_FILE = ROOT / "Data" / "processed" / "uniprot_clean.csv"


def main():

    if not INPUT_FILE.exists():
        raise FileNotFoundError(f"UniProt file not found: {INPUT_FILE}")

    # load UniProt data
    uniprot = pd.read_csv(INPUT_FILE, sep="\t")

    # clean the protein sequences so formatting is consistent
    uniprot["Sequence"] = (
        uniprot["Sequence"]
        .astype(str)
        .str.strip()
        .str.upper()
    )

    # quick check that calculated sequence lengths match the UniProt length column
    if "Length" in uniprot.columns:
        seq_length = uniprot["Sequence"].str.len()
        reported_length = pd.to_numeric(uniprot["Length"], errors="coerce")

        mismatches = uniprot[seq_length != reported_length]

        print(f"rows: {len(uniprot)}")
        print(f"sequence length mismatches: {len(mismatches)}")

        if not mismatches.empty:
            print(
                mismatches[
                    ["Entry", "Entry Name", "Length", "Sequence"]
                ].head()
            )

    # save cleaned version for the later mapping scripts
    OUTPUT_FILE.parent.mkdir(parents=True, exist_ok=True)
    uniprot.to_csv(OUTPUT_FILE, index=False)


if __name__ == "__main__":
    main()