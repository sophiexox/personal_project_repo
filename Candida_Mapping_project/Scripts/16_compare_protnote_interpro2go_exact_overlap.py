from pathlib import Path
import re

import pandas as pd

# file paths
ROOT = Path(__file__).resolve().parents[1]

PROTNOTE_FILE = ROOT / "Data" / "raw" / "protnote" / "protnote_interpro2go_benchmark_predictions.tsv"
BENCHMARK_FILE = ROOT / "Results" / "14_interpro2go_benchmark_proteins.csv"

OUTPUT_SUMMARY = ROOT / "Results" / "16_protnote_interpro2go_exact_overlap_summary.csv"
OUTPUT_PROTEIN = ROOT / "Results" / "16_protnote_interpro2go_exact_overlap_by_protein.csv"
OUTPUT_MATCHES = ROOT / "Results" / "16_protnote_interpro2go_exact_matching_terms.csv"

THRESHOLDS = [0.5, 0.7, 0.8, 0.9]

GO_COLUMNS = [
    "Computed GO Component IDs",
    "Computed GO Function IDs",
    "Computed GO Process IDs",
]

GO_PATTERN = re.compile(r"GO:\d{7}")


def extract_go_terms(value):
    # return all GO accessions found in one FungiDB field
    if pd.isna(value):
        return set()

    return set(GO_PATTERN.findall(str(value)))


def combine_interpro_terms(row):
    terms = set()

    for column in GO_COLUMNS:
        terms.update(extract_go_terms(row[column]))

    return terms


def safe_ratio(numerator, denominator):
    return numerator / denominator if denominator else 0


def main():
    for path in [PROTNOTE_FILE, BENCHMARK_FILE]:
        if not path.exists():
            raise FileNotFoundError(path)

    protnote = pd.read_csv(PROTNOTE_FILE, sep="\t")
    benchmark = pd.read_csv(BENCHMARK_FILE, dtype=str)

    required = {"sequence_id", "rank", "go_id", "probability", "description"}
    missing = required - set(protnote.columns)

    if missing:
        raise KeyError(f"Missing ProtNote columns: {', '.join(sorted(missing))}")

    missing_go = [column for column in GO_COLUMNS if column not in benchmark.columns]
    if missing_go:
        raise KeyError(f"Missing benchmark GO columns: {', '.join(missing_go)}")

    protnote["probability"] = pd.to_numeric(
        protnote["probability"],
        errors="coerce",
    )

    if protnote["probability"].isna().any():
        raise ValueError("Non-numeric ProtNote probabilities found")

    # build one InterPro2GO term set for each benchmark protein
    benchmark["InterPro2GO_terms"] = benchmark.apply(
        combine_interpro_terms,
        axis=1,
    )

    interpro_lookup = dict(
        zip(benchmark["Gene ID"], benchmark["InterPro2GO_terms"])
    )

    protnote_ids = set(protnote["sequence_id"].dropna().astype(str))
    benchmark_ids = set(benchmark["Gene ID"].dropna().astype(str))
    common_ids = protnote_ids & benchmark_ids

    summary_rows = []
    protein_rows = []
    match_rows = []

    for threshold in THRESHOLDS:
        filtered = protnote[
            (protnote["probability"] >= threshold)
            & protnote["sequence_id"].isin(common_ids)
        ].copy()

        # make one ProtNote GO set per protein at this threshold
        protnote_sets = (
            filtered
            .groupby("sequence_id")["go_id"]
            .apply(set)
            .to_dict()
        )

        proteins_with_predictions = len(protnote_sets)
        proteins_with_match = 0
        proteins_without_match = 0
        total_protnote_terms = 0
        total_interpro_terms = 0
        total_matches = 0

        for gene_id in sorted(common_ids):
            pn_terms = protnote_sets.get(gene_id, set())
            ip_terms = interpro_lookup.get(gene_id, set())
            matches = pn_terms & ip_terms

            pn_count = len(pn_terms)
            ip_count = len(ip_terms)
            match_count = len(matches)

            if pn_count:
                total_protnote_terms += pn_count
                total_interpro_terms += ip_count

                if matches:
                    proteins_with_match += 1
                else:
                    proteins_without_match += 1

            protein_rows.append({
                "Gene ID": gene_id,
                "Threshold": threshold,
                "ProtNote Term Count": pn_count,
                "InterPro2GO Term Count": ip_count,
                "Exact Match Count": match_count,
                "Has ProtNote Prediction": pn_count > 0,
                "Has Exact Match": match_count > 0,
                "Exact Precision": safe_ratio(match_count, pn_count) if pn_count else None,
                "Exact Recall": safe_ratio(match_count, ip_count) if ip_count else None,
                "ProtNote GO Terms": ";".join(sorted(pn_terms)),
                "InterPro2GO GO Terms": ";".join(sorted(ip_terms)),
                "Exact Matching GO Terms": ";".join(sorted(matches)),
            })

            if matches:
                matching_predictions = filtered[
                    (filtered["sequence_id"] == gene_id)
                    & filtered["go_id"].isin(matches)
                ]

                for prediction in matching_predictions.itertuples(index=False):
                    match_rows.append({
                        "Gene ID": gene_id,
                        "Threshold": threshold,
                        "GO ID": prediction.go_id,
                        "ProtNote Probability": prediction.probability,
                        "ProtNote Rank": prediction.rank,
                        "Description": prediction.description,
                    })

            total_matches += match_count

        summary_rows.append({
            "Threshold": threshold,
            "Benchmark Proteins": len(common_ids),
            "Proteins With ProtNote Predictions": proteins_with_predictions,
            "Proteins With >=1 Exact Match": proteins_with_match,
            "Proteins With Predictions But No Exact Match": proteins_without_match,
            "Protein Exact Agreement Rate": safe_ratio(
                proteins_with_match,
                proteins_with_predictions,
            ),
            "ProtNote GO Terms Predicted": total_protnote_terms,
            "InterPro2GO Terms For Predicted Proteins": total_interpro_terms,
            "Exact Matching GO Terms": total_matches,
            "Term-Level Exact Precision": safe_ratio(
                total_matches,
                total_protnote_terms,
            ),
            "Term-Level Exact Recall": safe_ratio(
                total_matches,
                total_interpro_terms,
            ),
        })

    summary = pd.DataFrame(summary_rows)
    protein_results = pd.DataFrame(protein_rows)
    matches = pd.DataFrame(match_rows)

    summary.to_csv(OUTPUT_SUMMARY, index=False)
    protein_results.to_csv(OUTPUT_PROTEIN, index=False)
    matches.to_csv(OUTPUT_MATCHES, index=False)

    print(f"benchmark proteins with ProtNote output: {len(common_ids)}")
    print(f"benchmark proteins without ProtNote output: {len(benchmark_ids - protnote_ids)}")
    print(f"ProtNote proteins absent from benchmark: {len(protnote_ids - benchmark_ids)}")

    for result in summary_rows:
        print(
            f"threshold >= {result['Threshold']:.1f}: "
            f"{result['Proteins With ProtNote Predictions']} proteins predicted, "
            f"{result['Proteins With >=1 Exact Match']} with exact match, "
            f"agreement {result['Protein Exact Agreement Rate']:.1%}"
        )


if __name__ == "__main__":
    main()


# GAI declaration
# GAI was used around lines 99-116 when grouping ProtNote predictions into one
# set of GO terms per protein and comparing this against another GO-term set.
# i asked for an example of exact overlap between grouped identifiers, and was
# shown:
#
# grouped = df.groupby("your_id")["your_term"].apply(set).to_dict()
# predicted = grouped.get(identifier, set())
# reference = reference_lookup.get(identifier, set())
# matches = predicted & reference
#
# i adapted this to ProtNote and InterPro2GO GO accessions at each threshold.

# lines 164-185 were also supported by GAI. i asked how to calculate several
# summary rates safely when a denominator might be zero. an example used a small
# helper such as `return a / b if b else 0`, which became safe_ratio() and is used
# for protein agreement, term-level precision and term-level recall.