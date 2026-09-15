from pathlib import Path

import pandas as pd


INPUT_FILE = Path("GenesByTaxon_Summary (1).txt")
OUTPUT_FILE = Path("cauris_annotation_clean.csv")
HYPOTHETICAL_WITH_GO_FILE = Path("hypothetical_with_go.csv")
HYPOTHETICAL_WITHOUT_GO_FILE = Path("hypothetical_without_go.csv")
FUNCTIONALLY_NAMED_FILE = Path("functionally_named_proteins.csv")
PRIORITISED_HYPOTHETICAL_FILE = Path("prioritised_hypothetical_with_go.csv")


def main() -> None:
    """Extract the annotation columns needed for AI-tool validation."""

    if not INPUT_FILE.exists():
        raise FileNotFoundError(
            f"Could not find {INPUT_FILE}. "
            "Place the input file in the same folder as this script."
        )

    # The uploaded file is tab-separated.
    annotations = pd.read_csv(
        INPUT_FILE,
        sep="\t",
        dtype=str,
        keep_default_na=False
    )

    wanted_columns = [
        "Gene ID",
        "source_id",
        "Product Description",
        "Computed GO Component IDs",
        "Computed GO Components",
        "Computed GO Function IDs",
        "Computed GO Functions",
        "Computed GO Process IDs",
        "Computed GO Processes",
        "Curated GO Component IDs",
        "Curated GO Components",
        "Curated GO Function IDs",
        "Curated GO Functions",
        "Curated GO Process IDs",
        "Curated GO Processes",
        "EC numbers",
        "EC numbers from OrthoMCL"
    ]

    missing_columns = [
        column for column in wanted_columns
        if column not in annotations.columns
    ]

    if missing_columns:
        raise ValueError(
            "The following expected columns were not found:\n"
            + "\n".join(missing_columns)
        )

    clean_table = annotations[wanted_columns].copy()

    # Add columns ready for later ProtNLM and ProtNote results.
    clean_table["ProtNLM prediction"] = ""
    clean_table["ProtNLM score"] = ""
    clean_table["ProtNote GO predictions"] = ""
    clean_table["GO2SUM summary"] = ""
    clean_table["Validation category"] = ""
    clean_table["Notes"] = ""

    clean_table.to_csv(OUTPUT_FILE, index=False)

    is_hypothetical = (
        clean_table["Product Description"]
        .str.strip()
        .str.lower()
        .eq("hypothetical protein")
    )

    has_computed_go_mask = (
        clean_table[
            [
                "Computed GO Component IDs",
                "Computed GO Function IDs",
                "Computed GO Process IDs"
            ]
        ]
        .apply(
            lambda column: ~column.str.strip().isin(["", "N/A"]),
            axis=0
        )
        .any(axis=1)
    )

    hypothetical_with_go = clean_table[is_hypothetical & has_computed_go_mask].copy()
    hypothetical_without_go = clean_table[is_hypothetical & ~has_computed_go_mask].copy()
    functionally_named = clean_table[~is_hypothetical].copy()

    def count_go_terms(value: str) -> int:
        """Count semicolon-separated GO identifiers, ignoring empty and N/A values."""
        cleaned = value.strip()
        if cleaned in {"", "N/A"}:
            return 0
        return len([term for term in cleaned.split(";") if term.strip()])

    hypothetical_with_go["GO component count"] = (
        hypothetical_with_go["Computed GO Component IDs"].apply(count_go_terms)
    )
    hypothetical_with_go["GO function count"] = (
        hypothetical_with_go["Computed GO Function IDs"].apply(count_go_terms)
    )
    hypothetical_with_go["GO process count"] = (
        hypothetical_with_go["Computed GO Process IDs"].apply(count_go_terms)
    )
    hypothetical_with_go["Total computed GO count"] = (
        hypothetical_with_go["GO component count"]
        + hypothetical_with_go["GO function count"]
        + hypothetical_with_go["GO process count"]
    )

    prioritised_hypothetical = hypothetical_with_go.sort_values(
        by=[
            "Total computed GO count",
            "GO function count",
            "GO process count",
            "Gene ID"
        ],
        ascending=[False, False, False, True]
    ).copy()

    hypothetical_with_go.to_csv(HYPOTHETICAL_WITH_GO_FILE, index=False)
    hypothetical_without_go.to_csv(HYPOTHETICAL_WITHOUT_GO_FILE, index=False)
    functionally_named.to_csv(FUNCTIONALLY_NAMED_FILE, index=False)
    prioritised_hypothetical.to_csv(PRIORITISED_HYPOTHETICAL_FILE, index=False)

    hypothetical_count = int(is_hypothetical.sum())
    has_computed_go = int(has_computed_go_mask.sum())

    print(f"Total proteins: {len(clean_table)}")
    print(f"Hypothetical proteins: {hypothetical_count}")
    print(f"Proteins with at least one computed GO annotation: {has_computed_go}")
    print(f"Hypothetical proteins with computed GO: {len(hypothetical_with_go)}")
    print(f"Hypothetical proteins without computed GO: {len(hypothetical_without_go)}")
    print(f"Functionally named proteins: {len(functionally_named)}")
    print(f"Clean table saved to: {OUTPUT_FILE.resolve()}")
    print(f"Hypothetical with GO saved to: {HYPOTHETICAL_WITH_GO_FILE.resolve()}")
    print(f"Hypothetical without GO saved to: {HYPOTHETICAL_WITHOUT_GO_FILE.resolve()}")
    print(f"Functionally named proteins saved to: {FUNCTIONALLY_NAMED_FILE.resolve()}")

    print(
        "Prioritised hypothetical proteins saved to: "
        f"{PRIORITISED_HYPOTHETICAL_FILE.resolve()}"
    )

    print("\nTop 10 hypothetical proteins by computed GO evidence:")
    print(
        prioritised_hypothetical[
            [
                "Gene ID",
                "GO component count",
                "GO function count",
                "GO process count",
                "Total computed GO count"
            ]
        ]
        .head(10)
        .to_string(index=False)
    )


if __name__ == "__main__":
    main()