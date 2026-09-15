from pathlib import Path

import pandas as pd

# file paths
ROOT = Path(__file__).resolve().parents[1]

CANDIDATE_FILE = ROOT / "Results" / "fungidb_rel70_annotation_poor_candidate_list.csv"
FULL_FASTA = ROOT / "Data" / "raw" / "cgd" / "C_auris_B8441_current_orf_trans_all.fasta"
OUTPUT_FASTA = ROOT / "Results" / "annotation_poor_526_proteins.fasta"
UNMATCHED_FILE = ROOT / "Results" / "11_annotation_poor_ids_not_found_in_fasta.csv"


def read_fasta(path):
    # store each FASTA record by its first header field, which is the B9J08 gene ID
    records = {}
    gene_id = None
    header = None
    sequence = []

    def save_record():
        if gene_id is not None:
            records[gene_id] = {
                "header": header,
                "sequence": "".join(sequence).upper(),
            }

    with path.open(encoding="utf-8") as fasta:
        for raw_line in fasta:
            line = raw_line.strip()

            if not line:
                continue

            if line.startswith(">"):
                save_record()
                header = line
                gene_id = line[1:].split()[0]
                sequence = []
            else:
                sequence.append(line)

    save_record()
    return records


def write_fasta(records, gene_ids, path):
    # write selected proteins and wrap sequence lines at 60 characters
    with path.open("w", encoding="utf-8") as output:
        for gene_id in gene_ids:
            record = records[gene_id]
            output.write(f"{record['header']}\n")

            sequence = record["sequence"]
            for start in range(0, len(sequence), 60):
                output.write(sequence[start:start + 60] + "\n")


def main():
    for path in [CANDIDATE_FILE, FULL_FASTA]:
        if not path.exists():
            raise FileNotFoundError(path)

    candidates = pd.read_csv(CANDIDATE_FILE, dtype=str)

    if "Gene ID" not in candidates.columns:
        raise KeyError("Gene ID column not found in candidate file")

    candidate_ids = set(
        candidates["Gene ID"]
        .dropna()
        .astype(str)
        .str.strip()
    )
    candidate_ids.discard("")

    records = read_fasta(FULL_FASTA)

    matched = sorted(candidate_ids & records.keys())
    unmatched = sorted(candidate_ids - records.keys())

    # check for any matched records with no sequence
    empty_sequences = [
        gene_id
        for gene_id in matched
        if not records[gene_id]["sequence"]
    ]

    OUTPUT_FASTA.parent.mkdir(exist_ok=True)
    write_fasta(records, matched, OUTPUT_FASTA)

    pd.DataFrame({"Gene ID": unmatched}).to_csv(
        UNMATCHED_FILE,
        index=False,
    )

    lengths = [len(records[gene_id]["sequence"]) for gene_id in matched]

    print(f"candidate IDs: {len(candidate_ids)}")
    print(f"matched FASTA records: {len(matched)}")
    print(f"unmatched IDs: {len(unmatched)}")
    print(f"empty matched sequences: {len(empty_sequences)}")

    if lengths:
        print(f"shortest protein: {min(lengths)} aa")
        print(f"longest protein: {max(lengths)} aa")
        print(f"mean protein length: {sum(lengths) / len(lengths):.1f} aa")


if __name__ == "__main__":
    main()


# GAI declaration
# GAI was used for lines 14-43 when working out how to parse a multi-line FASTA
# file into one record per protein. i asked how to keep the current header and
# collect sequence lines until the next ">" header appeared. the response showed
# a generic example using variables such as current_header and current_sequence
# together with a small save_record() function. i adapted this so each record is
# stored using its B9J08 gene ID as the dictionary key.

# i used GAI again for lines 47-56 when writing the output FASTA. my question was
# how to wrap an amino-acid sequence to 60 characters per line instead of writing
# the full sequence on one line. the example used
# range(0, len(your_sequence), 60) and slices such as
# your_sequence[start:start + 60], which i used in write_fasta().

# lines 79-80 were based on a GAI suggestion for comparing requested IDs against
# the identifiers available in a dictionary. the generic response used set
# intersection to return matches and set subtraction to return missing IDs.
# i applied this to candidate_ids and records.keys() to create the matched and
# unmatched candidate lists.