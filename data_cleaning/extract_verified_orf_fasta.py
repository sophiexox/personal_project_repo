from pathlib import Path


INPUT_FASTA = Path("C_auris_B8441_current_orf_trans_all.fasta")
OUTPUT_FASTA = Path("cauris_verified_orf_pilot_94.fasta")
HEADER_MARKER = "Verified ORF;"


def read_fasta(path: Path):
    """Yield FASTA records as (header, sequence)."""
    header = None
    sequence_lines = []

    with path.open("r", encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.strip()

            if not line:
                continue

            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(sequence_lines)

                header = line[1:]
                sequence_lines = []
            else:
                sequence_lines.append(line)

        if header is not None:
            yield header, "".join(sequence_lines)


def write_wrapped_sequence(handle, sequence: str, width: int = 60) -> None:
    """Write a protein sequence over fixed-width FASTA lines."""
    for start in range(0, len(sequence), width):
        handle.write(sequence[start:start + width] + "\n")


def main() -> None:
    if not INPUT_FASTA.exists():
        raise FileNotFoundError(
            f"Could not find {INPUT_FASTA}. "
            "Place the FASTA file in the same folder as this script."
        )

    total_records = 0
    verified_records = 0
    verified_ids = []

    with OUTPUT_FASTA.open("w", encoding="utf-8") as output_handle:
        for header, sequence in read_fasta(INPUT_FASTA):
            total_records += 1

            if HEADER_MARKER not in header:
                continue

            verified_records += 1
            gene_id = header.split()[0]
            verified_ids.append(gene_id)

            output_handle.write(f">{header}\n")
            write_wrapped_sequence(output_handle, sequence)

    print(f"Total FASTA records examined: {total_records}")
    print(f"Verified ORF records extracted: {verified_records}")
    print(f"Output FASTA: {OUTPUT_FASTA.resolve()}")

    if verified_records != 94:
        print(
            f"Warning: expected 94 Verified ORF records, "
            f"but extracted {verified_records}."
        )

    if len(verified_ids) != len(set(verified_ids)):
        print("Warning: duplicate gene IDs were detected.")
    else:
        print("All extracted gene IDs are unique.")


if __name__ == "__main__":
    main()