from pathlib import Path
import re

import pandas as pd

# file paths
ROOT = Path(__file__).resolve().parents[1]

PROTNOTE_FILE = ROOT / "Data" / "raw" / "protnote" / "protnote_interpro2go_benchmark_predictions.tsv"
BENCHMARK_FILE = ROOT / "Results" / "14_interpro2go_benchmark_proteins.csv"

OUTPUT_SUMMARY = ROOT / "Results" / "17_protnote_interpro2go_overlap_by_namespace.csv"
OUTPUT_PROTEIN = ROOT / "Results" / "17_protnote_interpro2go_overlap_by_namespace_by_protein.csv"

THRESHOLDS = [0.5, 0.7, 0.8, 0.9]

NAMESPACES = {
    "CC": "Computed GO Component IDs",
    "MF": "Computed GO Function IDs",
    "BP": "Computed GO Process IDs",
}

GO_PATTERN = re.compile(r"GO:\d{7}")


def extract_go_terms(value):
    # return all GO accessions found in one annotation field
    if pd.isna(value):
        return set()

    return set(GO_PATTERN.findall(str(value)))


def safe_ratio(numerator, denominator):
    return numerator / denominator if denominator else 0


def main():
    for path in [PROTNOTE_FILE, BENCHMARK_FILE]:
        if not path.exists():
            raise FileNotFoundError(path)

    protnote = pd.read_csv(PROTNOTE_FILE, sep="\t")
    benchmark = pd.read_csv(BENCHMARK_FILE, dtype=str)

    needed = {"sequence_id", "go_id", "probability"}
    missing = needed - set(protnote.columns)

    if missing:
        raise KeyError(f"Missing ProtNote columns: {', '.join(sorted(missing))}")

    missing_go = [
        column for column in NAMESPACES.values()
        if column not in benchmark.columns
    ]
    if missing_go:
        raise KeyError(f"Missing benchmark GO columns: {', '.join(missing_go)}")

    protnote["probability"] = pd.to_numeric(
        protnote["probability"],
        errors="raise",
    )

    # build namespace-specific InterPro2GO term sets
    for namespace, column in NAMESPACES.items():
        benchmark[f"{namespace}_terms"] = benchmark[column].apply(extract_go_terms)

    # infer GO namespace from the benchmark annotations
    go_to_namespace = {}
    ambiguous = set()

    for namespace in NAMESPACES:
        for terms in benchmark[f"{namespace}_terms"]:
            for go_id in terms:
                if go_id in go_to_namespace and go_to_namespace[go_id] != namespace:
                    ambiguous.add(go_id)
                else:
                    go_to_namespace[go_id] = namespace

    for go_id in ambiguous:
        go_to_namespace.pop(go_id, None)

    protnote["namespace"] = protnote["go_id"].map(go_to_namespace)

    protnote_ids = set(protnote["sequence_id"].dropna().astype(str))
    benchmark_ids = set(benchmark["Gene ID"].dropna().astype(str))
    common_ids = protnote_ids & benchmark_ids

    interpro_lookup = {
        namespace: dict(zip(benchmark["Gene ID"], benchmark[f"{namespace}_terms"]))
        for namespace in NAMESPACES
    }

    summary_rows = []
    protein_rows = []

    for threshold in THRESHOLDS:
        threshold_data = protnote[
            (protnote["probability"] >= threshold)
            & protnote["sequence_id"].isin(common_ids)
        ]

        for namespace in NAMESPACES:
            ns_predictions = threshold_data[
                threshold_data["namespace"] == namespace
            ]

            pn_lookup = (
                ns_predictions
                .groupby("sequence_id")["go_id"]
                .apply(set)
                .to_dict()
            )

            # only proteins with reference terms in this namespace are eligible
            eligible_ids = {
                gene_id
                for gene_id in common_ids
                if interpro_lookup[namespace].get(gene_id, set())
            }

            proteins_with_prediction = 0
            proteins_with_match = 0
            total_pn_terms = 0
            total_ip_terms = 0
            total_matches = 0

            for gene_id in sorted(eligible_ids):
                pn_terms = pn_lookup.get(gene_id, set())
                ip_terms = interpro_lookup[namespace].get(gene_id, set())
                matches = pn_terms & ip_terms

                if pn_terms:
                    proteins_with_prediction += 1
                    total_pn_terms += len(pn_terms)
                    total_ip_terms += len(ip_terms)
                    total_matches += len(matches)

                    if matches:
                        proteins_with_match += 1

                protein_rows.append({
                    "Gene ID": gene_id,
                    "Threshold": threshold,
                    "Namespace": namespace,
                    "ProtNote Term Count": len(pn_terms),
                    "InterPro2GO Term Count": len(ip_terms),
                    "Exact Match Count": len(matches),
                    "Has ProtNote Prediction": bool(pn_terms),
                    "Has Exact Match": bool(matches),
                    "ProtNote GO Terms": ";".join(sorted(pn_terms)),
                    "InterPro2GO GO Terms": ";".join(sorted(ip_terms)),
                    "Exact Matching GO Terms": ";".join(sorted(matches)),
                })

            summary_rows.append({
                "Threshold": threshold,
                "Namespace": namespace,
                "Eligible Proteins": len(eligible_ids),
                "Proteins With ProtNote Predictions": proteins_with_prediction,
                "Proteins With >=1 Exact Match": proteins_with_match,
                "Protein Exact Agreement Rate": safe_ratio(
                    proteins_with_match,
                    proteins_with_prediction,
                ),
                "ProtNote GO Terms": total_pn_terms,
                "InterPro2GO Terms": total_ip_terms,
                "Exact Matching GO Terms": total_matches,
                "Term-Level Exact Precision": safe_ratio(
                    total_matches,
                    total_pn_terms,
                ),
                "Term-Level Exact Recall": safe_ratio(
                    total_matches,
                    total_ip_terms,
                ),
            })

    summary = pd.DataFrame(summary_rows)
    protein_results = pd.DataFrame(protein_rows)

    summary.to_csv(OUTPUT_SUMMARY, index=False)
    protein_results.to_csv(OUTPUT_PROTEIN, index=False)

    print(f"proteins present in both datasets: {len(common_ids)}")
    print(f"ProtNote rows mapped to namespace: {protnote['namespace'].notna().sum()}")
    print(f"ProtNote rows without namespace: {protnote['namespace'].isna().sum()}")
    print(f"ambiguous GO IDs excluded from namespace lookup: {len(ambiguous)}")

    for row in summary_rows:
        print(
            f"{row['Namespace']} >= {row['Threshold']:.1f}: "
            f"{row['Proteins With ProtNote Predictions']} predicted, "
            f"{row['Proteins With >=1 Exact Match']} matched, "
            f"agreement {row['Protein Exact Agreement Rate']:.1%}"
        )


if __name__ == "__main__":
    main()


# GAI declaration
# GAI was used for lines 68-83 when building a GO ID to namespace lookup.
# i asked how to assign values from several groups to a dictionary while catching
# identifiers that appeared in more than one group. an example used:
#
# if item in lookup and lookup[item] != group:
#     ambiguous.add(item)
# else:
#     lookup[item] = group
#
# i adapted this to identify CC, MF and BP terms and remove ambiguous GO IDs.

# for lines 103-120, i asked how to group prediction rows into a set of terms per
# protein and then restrict analysis to IDs that had reference terms available.
# the example used:
#
# lookup = df.groupby("your_id")["your_term"].apply(set).to_dict()
# eligible = {x for x in ids if reference_lookup.get(x, set())}
#
# this was adapted for each GO namespace and ProtNote probability threshold.