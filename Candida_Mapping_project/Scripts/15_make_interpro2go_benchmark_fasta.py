from pathlib import Path
from collections import defaultdict

import pandas as pd

# file paths
ROOT = Path(__file__).resolve().parents[1]

BENCHMARK_FILE = ROOT / "Results" / "14_interpro2go_benchmark_proteins.csv"
FULL_FASTA = ROOT / "Data" / "raw" / "cgd" / "C_auris_B8441_current_orf_trans_all.fasta"
OUTPUT_FASTA = ROOT / "Results" / "15_interpro2go_benchmark_proteins_clean.fasta"
UNMATCHED_FILE = ROOT / "Results" / "15_interpro2go_ids_not_found_in_fasta.csv"
DUPLICATE_FILE = ROOT / "Results" / "15_duplicate_benchmark_fasta_ids.csv"
NONSTANDARD_FILE = ROOT / "Results" / "15_nonstandard_residue_exclusions.csv"
CONFLICT_FILE = ROOT / "Results" / "15_conflicting_duplicate_sequences.csv"

ALLOWED_AA = set("ACDEFGHIKLMNPQRSTVWY")


def read_fasta(path, wanted_ids):
    # keep every sequence found for benchmark IDs so duplicates can be checked
    records = defaultdict(list)
    gene_id = None
    sequence = []

    def save_record():
        if gene_id in wanted_ids:
            records[gene_id].append("".join(sequence).upper())

    with path.open(encoding="utf-8") as fasta:
        for raw_line in fasta:
            line = raw_line.strip()

            if not line:
                continue

            if line.startswith(">"):
                save_record()
                gene_id = line[1:].split()[0]
                sequence = []
            else:
                sequence.append(line)

    save_record()
    return records


def check_duplicates(records):
    # separate safe duplicate records from conflicting sequences
    duplicate_rows = []
    conflict_rows = []
    unique_records = {}

    for gene_id, sequences in records.items():
        unique_sequences = list(dict.fromkeys(sequences))

        if len(sequences) > 1:
            duplicate_rows.append({
                "Gene ID": gene_id,
                "FASTA record count": len(sequences),
                "Unique sequence count": len(unique_sequences),
                "Sequences identical": len(unique_sequences) == 1,
            })

        if len(unique_sequences) == 1:
            unique_records[gene_id] = unique_sequences[0]
        else:
            conflict_rows.append({
                "Gene ID": gene_id,
                "FASTA record count": len(sequences),
                "Unique sequence count": len(unique_sequences),
                "Sequence lengths": ";".join(
                    str(len(sequence))
                    for sequence in sorted(unique_sequences, key=len)
                ),
            })

    duplicates = pd.DataFrame(duplicate_rows)
    conflicts = pd.DataFrame(conflict_rows)

    return unique_records, duplicates, conflicts


def remove_nonstandard(records):
    # ProtNote input is restricted to the 20 standard amino acids
    clean = {}
    excluded = []

    for gene_id, sequence in records.items():
        unsupported = sorted(set(sequence) - ALLOWED_AA)

        if unsupported:
            excluded.append({
                "Gene ID": gene_id,
                "Unsupported residues": ",".join(unsupported),
                "Protein length": len(sequence),
            })
        else:
            clean[gene_id] = sequence

    return clean, pd.DataFrame(excluded)


def write_fasta(records, path):
    with path.open("w", encoding="utf-8") as output:
        for gene_id in sorted(records):
            sequence = records[gene_id]
            output.write(f">{gene_id}\n")

            for start in range(0, len(sequence), 60):
                output.write(sequence[start:start + 60] + "\n")


def main():
    for path in [BENCHMARK_FILE, FULL_FASTA]:
        if not path.exists():
            raise FileNotFoundError(path)

    benchmark = pd.read_csv(BENCHMARK_FILE, dtype=str)

    if "Gene ID" not in benchmark.columns:
        raise KeyError("Gene ID column not found in benchmark file")

    benchmark_ids = set(
        benchmark["Gene ID"]
        .dropna()
        .astype(str)
        .str.strip()
    )
    benchmark_ids.discard("")

    records = read_fasta(FULL_FASTA, benchmark_ids)

    matched_ids = set(records)
    unmatched_ids = sorted(benchmark_ids - matched_ids)

    unique_records, duplicates, conflicts = check_duplicates(records)
    clean_records, nonstandard = remove_nonstandard(unique_records)

    OUTPUT_FASTA.parent.mkdir(exist_ok=True)
    write_fasta(clean_records, OUTPUT_FASTA)

    pd.DataFrame({"Gene ID": unmatched_ids}).to_csv(
        UNMATCHED_FILE,
        index=False,
    )
    duplicates.to_csv(DUPLICATE_FILE, index=False)
    nonstandard.to_csv(NONSTANDARD_FILE, index=False)
    conflicts.to_csv(CONFLICT_FILE, index=False)

    lengths = [len(sequence) for sequence in clean_records.values()]

    print(f"benchmark proteins: {len(benchmark_ids)}")
    print(f"matched IDs: {len(matched_ids)}")
    print(f"unmatched IDs: {len(unmatched_ids)}")
    print(f"duplicate FASTA IDs: {len(duplicates)}")
    print(f"conflicting duplicate sequences: {len(conflicts)}")
    print(f"non-standard residue exclusions: {len(nonstandard)}")
    print(f"final ProtNote-compatible proteins: {len(clean_records)}")

    if lengths:
        print(f"shortest protein: {min(lengths)} aa")
        print(f"longest protein: {max(lengths)} aa")
        print(f"mean protein length: {sum(lengths) / len(lengths):.1f} aa")


if __name__ == "__main__":
    main()


# GAI declaration
# GAI was used for lines 20-45 to help with reading the FASTA while retaining
# every sequence associated with the same gene ID. i needed duplicate IDs to be
# kept rather than overwritten, so i asked how to parse a FASTA into a dictionary
# where each ID can contain multiple sequences. the example response was similar to:
#
# from collections import defaultdict
#
# records = defaultdict(list)
# current_id = None
# sequence_parts = []
#
# def save_record():
#     if current_id is not None:
#         sequence = "".join(sequence_parts)
#         records[current_id].append(sequence)
#
# with open("your_file.fasta") as fasta:
#     for line in fasta:
#         line = line.strip()
#
#         if line.startswith(">"):
#             save_record()
#             current_id = line[1:].split()[0]
#             sequence_parts = []
#         else:
#             sequence_parts.append(line)
#
# save_record()
#
# i adapted this so only IDs present in benchmark_ids are stored, and sequences
# are converted to uppercase before later QC.

# lines 48-82 also used GAI assistance. i asked how to distinguish duplicate IDs
# where every sequence is identical from duplicates where the sequences actually
# disagree, without arbitrarily selecting one sequence. the response suggested
# reducing each list to its unique values and checking its length, for example:
#
# for identifier, sequences in records.items():
#     unique_sequences = list(dict.fromkeys(sequences))
#
#     if len(sequences) > 1:
#         print(identifier, "has duplicate records")
#
#     if len(unique_sequences) == 1:
#         clean_records[identifier] = unique_sequences[0]
#     else:
#         conflicting_ids.append(identifier)
#
# i expanded this to save the FASTA record count, number of unique sequences and
# sequence lengths. only IDs with one unique sequence are carried forward.

# GAI was used again for lines 85-102 when checking whether sequences contained
# amino-acid characters outside the standard 20 residues. i asked for a simple
# way to find characters present in a sequence but absent from an allowed set.
# the generated example used set subtraction:
#
# allowed = set("ACDEFGHIKLMNPQRSTVWY")
#
# for identifier, sequence in your_sequences.items():
#     unsupported = sorted(set(sequence) - allowed)
#
#     if unsupported:
#         print(identifier, unsupported)
#     else:
#         clean_sequences[identifier] = sequence
#
# i used this approach in remove_nonstandard() and recorded excluded IDs,
# unsupported residues and protein length in a separate QC table.

# for lines 105-112, i asked how to write sequences back to FASTA with a fixed
# line width rather than placing each complete protein sequence on one line.
# ChatGPT supplied a general example:
#
# with open("your_output.fasta", "w") as output:
#     for identifier, sequence in your_records.items():
#         output.write(f">{identifier}\\n")
#
#         for start in range(0, len(sequence), 60):
#             output.write(sequence[start:start + 60] + "\\n")
#
# this was adapted in write_fasta(), with IDs sorted first so the generated
# benchmark FASTA has a consistent order.

# lines 135-136 use a smaller GAI suggestion for checking FASTA coverage.
# i asked how to identify expected IDs that were missing after parsing the file.
# the response demonstrated set subtraction:
#
# matched_ids = set(records)
# missing_ids = expected_ids - matched_ids
#
# i applied this to benchmark_ids and the IDs recovered from records, then saved
# the missing values to 15_interpro2go_ids_not_found_in_fasta.csv.