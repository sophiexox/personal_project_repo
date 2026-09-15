from pathlib import Path
import re

import pandas as pd

# file paths
ROOT = Path(__file__).resolve().parents[1]

INPUT_FILE = ROOT / "All_auris_genes_beta_fungidb_rel70.tsv"
OUTPUT_FILE = ROOT / "Results" / "14_interpro2go_benchmark_proteins.csv"

GO_COLUMNS = {
    "Component": "Computed GO Component IDs",
    "Function": "Computed GO Function IDs",
    "Process": "Computed GO Process IDs",
}

GO_PATTERN = re.compile(r"GO:\d{7}")


def extract_go_terms(value):
    # pull GO accessions out of a FungiDB annotation field
    if pd.isna(value):
        return []

    return GO_PATTERN.findall(str(value))


def combine_unique(row):
    # combine the three namespaces without keeping duplicate GO IDs
    terms = []

    for column in [
        "InterPro2GO Component Terms",
        "InterPro2GO Function Terms",
        "InterPro2GO Process Terms",
    ]:
        for term in row[column]:
            if term not in terms:
                terms.append(term)

    return terms


def main():
    if not INPUT_FILE.exists():
        raise FileNotFoundError(INPUT_FILE)

    df = pd.read_csv(INPUT_FILE, sep="\t", dtype=str)

    needed = ["Gene ID", *GO_COLUMNS.values()]
    missing = [column for column in needed if column not in df.columns]

    if missing:
        raise KeyError(f"Missing columns: {', '.join(missing)}")

    # clean the three computed GO fields
    for column in GO_COLUMNS.values():
        df[column] = df[column].fillna("").str.strip()

    has_go = {
        namespace: df[column].ne("")
        for namespace, column in GO_COLUMNS.items()
    }

    any_go = pd.concat(has_go.values(), axis=1).any(axis=1)
    all_three = pd.concat(has_go.values(), axis=1).all(axis=1)

    # only proteins with at least one computed GO term enter the benchmark
    benchmark = df[any_go].copy()

    for namespace, column in GO_COLUMNS.items():
        benchmark[f"InterPro2GO {namespace} Terms"] = (
            benchmark[column].apply(extract_go_terms)
        )

    benchmark["InterPro2GO All Terms"] = benchmark.apply(
        combine_unique,
        axis=1,
    )
    benchmark["InterPro2GO Term Count"] = (
        benchmark["InterPro2GO All Terms"].apply(len)
    )

    # convert GO lists to text before writing the CSV
    term_columns = [
        "InterPro2GO Component Terms",
        "InterPro2GO Function Terms",
        "InterPro2GO Process Terms",
        "InterPro2GO All Terms",
    ]

    for column in term_columns:
        benchmark[column] = benchmark[column].apply(
            lambda terms: ";".join(terms)
        )

    output_columns = [
        "Gene ID",
        "source_id",
        "Product Description",
        "Protein Length",
        "PFam ID",
        "InterPro ID",
        *GO_COLUMNS.values(),
        *term_columns,
        "InterPro2GO Term Count",
    ]
    output_columns = [
        column for column in output_columns
        if column in benchmark.columns
    ]

    OUTPUT_FILE.parent.mkdir(exist_ok=True)
    benchmark[output_columns].to_csv(OUTPUT_FILE, index=False)

    unique_terms = {
        term
        for value in benchmark["InterPro2GO All Terms"]
        for term in value.split(";")
        if term
    }

    print(f"proteins with Component GO: {int(has_go['Component'].sum())}")
    print(f"proteins with Function GO: {int(has_go['Function'].sum())}")
    print(f"proteins with Process GO: {int(has_go['Process'].sum())}")
    print(f"proteins with any computed GO: {int(any_go.sum())}")
    print(f"proteins with all three namespaces: {int(all_three.sum())}")
    print(f"benchmark proteins: {len(benchmark)}")
    print(f"unique GO terms in benchmark: {len(unique_terms)}")


if __name__ == "__main__":
    main()


# GAI declaration
# GAI was used for lines 21-26 when extracting GO accessions from the FungiDB
# fields. i asked how to extract every identifier matching the format GO: followed
# by seven digits when the cell might also contain other text or delimiters.
# ChatGPT returned a generic example along the lines of:
#
# import re
#
# go_pattern = re.compile(r"GO:\d{7}")
#
# def extract_go_ids(value):
#     if value is None:
#         return []
#
#     value = str(value)
#     matches = go_pattern.findall(value)
#     return matches
#
# dataframe["GO IDs"] = dataframe["your_go_column"].apply(extract_go_ids)
#
# i used the same approach in extract_go_terms(), but added pd.isna() because the
# source data is being handled through pandas and may contain missing values.

# lines 29-42 were based on another GAI response. i needed to combine GO IDs from
# three list-valued columns, remove duplicates, but keep the order in which each
# term first appeared. the suggested generic code was:
#
# columns_to_combine = [
#     "your_column_1",
#     "your_column_2",
#     "your_column_3",
# ]
#
# def combine_unique(row):
#     combined = []
#
#     for column in columns_to_combine:
#         for item in row[column]:
#             if item not in combined:
#                 combined.append(item)
#
#     return combined
#
# dataframe["combined_terms"] = dataframe.apply(
#     combine_unique,
#     axis=1
# )
#
# i adapted the column names to the Component, Function and Process InterPro2GO
# lists and used the result as the combined InterPro2GO All Terms field.

# GAI was also used for lines 58-67. i asked how to make one mask for proteins
# that had a value in at least one of three GO columns, and another for proteins
# where all three columns contained a value. ChatGPT suggested something like:
#
# masks = {
#     "group_1": dataframe["your_column_1"].ne(""),
#     "group_2": dataframe["your_column_2"].ne(""),
#     "group_3": dataframe["your_column_3"].ne(""),
# }
#
# mask_table = pd.concat(masks.values(), axis=1)
#
# has_any = mask_table.any(axis=1)
# has_all = mask_table.all(axis=1)
#
# filtered = dataframe[has_any].copy()
#
# i adapted this to the Component, Function and Process GO namespaces. keeping the
# masks in a dictionary also lets the individual namespace coverage be reported.

# lines 117-122 use a set comprehension that came from a GAI query about counting
# distinct identifiers stored as semicolon-separated text across many dataframe
# rows. the generic solution given was:
#
# unique_ids = {
#     identifier
#     for cell in dataframe["your_combined_column"]
#     for identifier in str(cell).split(";")
#     if identifier
# }
#
# print(len(unique_ids))
#
# i used this on InterPro2GO All Terms to calculate the total number of unique GO
# accessions represented in the benchmark dataset.