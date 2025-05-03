# Ensure you have Biopython installed:
# pip install biopython
# vibe coded with gemini 2.5 pro

import os
import sys

from Bio import AlignIO

for file in os.listdir("/home/delalamo/sabr/a3ms/"):
    if not file.endswith("a3m"):
        continue
    input_a3m_file = os.path.join("/home/delalamo/sabr/a3ms/", file)
    output_sto_file = os.path.join("/home/delalamo/sabr/stos/", file[:-4] + ".sto")
    # --- Configuration ---
    # input_a3m_file = "/home/delalamo/sabr/a3ms/human_light.a3m"  # Replace with your A3M file path
    # output_sto_file = "/home/delalamo/sabr/stos/human_heavy.sto" # Replace with your desired output Stockholm file path
    # --- End Configuration ---

    print(f"Attempting to convert '{input_a3m_file}' (A3M) to '{output_sto_file}' (Stockholm)...")

    try:
        # --- Read the A3M file as FASTA ---
        # Biopython reads A3M as FASTA format. This generally works but loses
        # the distinction between uppercase (match) and lowercase (insert) states
        # specific to the A3M format definition used by tools like HHblits.
        print(f"Reading A3M file '{input_a3m_file}' using 'fasta' format...")

        # Use AlignIO.read if your A3M file contains only ONE alignment
        alignment = AlignIO.read(input_a3m_file, "fasta")
        print(f"Read {len(alignment)} sequences with alignment length {alignment.get_alignment_length()}.")

        # --- Alternatively: Use AlignIO.parse if your file might contain MULTIPLE alignments ---
        # alignments = list(AlignIO.parse(input_a3m_file, "fasta"))
        # if not alignments:
        #     raise ValueError("No alignments found in the input file.")
        # print(f"Read {len(alignments)} alignment(s). Processing the first one.")
        # alignment = alignments[0] # Or loop through 'alignments' if needed

        # --- Write the alignment object(s) to Stockholm format ---
        print(f"Writing alignment to Stockholm file '{output_sto_file}'...")

        # Use AlignIO.write for a single alignment object (from AlignIO.read)
        count = AlignIO.write(alignment, output_sto_file, "stockholm")

        # --- If using AlignIO.parse for multiple alignments: ---
        # count = AlignIO.write(alignments, output_sto_file, "stockholm")

        print(f"Successfully converted '{input_a3m_file}' to '{output_sto_file}'.")
        # The 'count' usually indicates the number of alignment blocks written (typically 1).
        # print(f"Wrote {count} alignment block(s).")

    except FileNotFoundError:
        print(f"Error: Input file '{input_a3m_file}' not found.", file=sys.stderr)
        sys.exit(1)
    except ValueError as e:
        # Catches errors like format issues during parsing or writing
        print(f"Error during conversion: {e}", file=sys.stderr)
        print("Please ensure the input file is a valid A3M/FASTA-like format and Biopython can parse it.", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        # Catch any other unexpected errors
        print(f"An unexpected error occurred: {e}", file=sys.stderr)
        sys.exit(1)