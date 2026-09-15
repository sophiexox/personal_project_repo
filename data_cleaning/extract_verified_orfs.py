# extract_verified_orfs.py

"""
Extract all FASTA records containing 'Verified ORF'
from the Candida auris B8441 protein FASTA.

Input:
    C_auris_B8441_current_orf_trans_all.fasta

Output:
    verified_orfs_only.fasta
"""

input_fasta = "C_auris_B8441_current_orf_trans_all.fasta"
output_fasta = "verified_orfs_only.fasta"

# Read FASTA file
with open(input_fasta, "r") as infile:
    lines = infile.readlines()

# Find header positions
header_positions = [
    i for i, line in enumerate(lines)
    if line.startswith(">")
]

verified_records = []
verified_count = 0

for idx, start in enumerate(header_positions):

    header = lines[start]

    # Find end of current FASTA record
    if idx < len(header_positions) - 1:
        end = header_positions[idx + 1]
    else:
        end = len(lines)

    # Keep only Verified ORFs
    if "Verified ORF" in header:

        verified_count += 1

        record = lines[start:end]
        verified_records.extend(record)

# Write output FASTA
with open(output_fasta, "w") as outfile:
    outfile.writelines(verified_records)

print(f"Verified ORFs extracted: {verified_count}")
print(f"Saved to: {output_fasta}")