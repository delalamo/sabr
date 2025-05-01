#!/usr/bin/env python3

import argparse
import os
from collections import defaultdict

# Biopython is required: pip install biopython
try:
    import warnings

    from Bio.PDB import PDBParser
    from Bio.PDB.PDBExceptions import PDBConstructionWarning
    # Removed three_to_one import
    from Bio.PDB.Polypeptide import is_aa

    # Suppress specific Biopython warnings if desired
    warnings.simplefilter('ignore', PDBConstructionWarning)
except ImportError:
    print("Error: Biopython library not found.")
    print("Please install it using: pip install biopython")
    exit(1)
except Exception as e:
    print(f"Error importing Biopython components: {e}")
    exit(1)


# --- Configuration ---
# Default input directory (containing the subdirectories like 'human_heavy')
DEFAULT_INPUT_BASE_DIR = './pdb_segments'
# Expected subdirectories to process
DEFAULT_SUBDIRS = [
    'human_heavy', 'human_light', 'mouse_heavy',
    'mouse_light', 'camelid_heavy'
]
# Target sequence length for the alignment
TARGET_SEQUENCE_LENGTH = 128
# Gap character
GAP_CHAR = '-'
# Unknown amino acid character
UNKNOWN_AA_CHAR = 'X'

# --- Amino Acid Code Mapping ---
# Manually define the standard three-to-one mapping
THREE_TO_ONE_MAP = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V"
}


# --- Helper Functions ---

def get_sequence_identifier(filename):
    """
    Extracts a sequence identifier (e.g., 'pdbcode_chainid') from the filename.
    Assumes filename format like 'pdbcode_chainid_start-end.pdb'.
    """
    base = os.path.basename(filename)
    parts = base.split('_')
    if len(parts) >= 2:
        return f"{parts[0]}_{parts[1]}"
    else:
        # Fallback if filename format is unexpected
        return os.path.splitext(base)[0]

def extract_aligned_sequence(pdb_filepath, max_length=TARGET_SEQUENCE_LENGTH):
    """
    Parses a PDB file, extracts the amino acid sequence from the first chain,
    aligns it to a fixed length based on residue number, ignoring insertion codes.

    Args:
        pdb_filepath (str): Path to the PDB file.
        max_length (int): The target length of the sequence (e.g., 128).

    Returns:
        tuple: (sequence_identifier, sequence_string) or (None, None) on error.
    """
    parser = PDBParser(QUIET=True)
    identifier = get_sequence_identifier(pdb_filepath)
    structure_id = os.path.splitext(os.path.basename(pdb_filepath))[0]

    try:
        structure = parser.get_structure(structure_id, pdb_filepath)
    except Exception as e:
        # More specific error catching could be added (e.g., FileNotFoundError)
        print(f"  Warning: Could not parse PDB file {pdb_filepath}: {e}")
        return None, None

    # Assume the first model and first chain (segment files should be simple)
    if not structure or len(structure) == 0:
        print(f"  Warning: No models found in {pdb_filepath}")
        return None, None
    model = structure[0]

    if not model or len(model.child_list) == 0:
         print(f"  Warning: No chains found in model 0 of {pdb_filepath}")
         return None, None
    # Ensure we get the chain object correctly
    # Segment files *should* only have one chain after the previous script
    try:
        chain = next(model.get_chains()) # More robust way to get the first chain
    except StopIteration:
         print(f"  Warning: No chains found using iterator in model 0 of {pdb_filepath}")
         return None, None


    # Store residues by their sequence number (resseq)
    residues_dict = {}
    for residue in chain:
        # residue.id is a tuple: (hetfield, resseq, icode)
        hetfield, resseq, icode = residue.id

        # Skip residues with insertion codes
        if icode.strip(): # Check if icode is not empty or whitespace
            continue

        # Skip non-standard residues or HETATMs based on hetfield
        # Also check if it's a standard amino acid using is_aa for robustness
        if hetfield.strip() or not is_aa(residue, standard=True):
             continue

        # Only consider residues within the target range (1 to max_length)
        if 1 <= resseq <= max_length:
            resname = residue.get_resname().upper() # Use upper case for map lookup
            # Get the one-letter code using the defined map
            aa_code = THREE_TO_ONE_MAP.get(resname)

            if aa_code:
                # Store only the first encountered residue for a given number
                # (although insertion codes are skipped, this is a safeguard)
                if resseq not in residues_dict:
                    residues_dict[resseq] = aa_code
            # else: # Optional: Handle non-standard 3-letter codes if needed
                # print(f"  Warning: Unknown residue name '{resname}' in {pdb_filepath} at position {resseq}. Skipping.")


    # Build the final sequence string of target length
    sequence = []
    for i in range(1, max_length + 1):
        sequence.append(residues_dict.get(i, GAP_CHAR)) # Use gap if residue number is missing

    return identifier, "".join(sequence)


# --- Main Processing Logic ---

def create_a3m_alignments(input_base_dir, subdirs_to_process):
    """
    Processes PDB segment files in specified subdirectories and creates A3M alignments.
    """
    print("--- PDB Segments to A3M Alignment Script ---")
    print(f"Input base directory: {input_base_dir}")
    print(f"Processing subdirectories: {', '.join(subdirs_to_process)}")
    print(f"Target sequence length: {TARGET_SEQUENCE_LENGTH}")

    if not os.path.isdir(input_base_dir):
        print(f"Error: Input base directory not found: {input_base_dir}")
        return

    # Process each specified subdirectory
    for subdir_name in subdirs_to_process:
        subdir_path = os.path.join(input_base_dir, subdir_name)
        # Output A3M file in the base directory
        a3m_output_file = os.path.join(input_base_dir, f"{subdir_name}.a3m")

        print(f"\nProcessing subdirectory: {subdir_path}...")

        if not os.path.isdir(subdir_path):
            print(f"  Warning: Subdirectory '{subdir_name}' not found. Skipping.")
            continue

        # List comprehension to find PDB files, case-insensitive
        pdb_files = [f for f in os.listdir(subdir_path)
                     if os.path.isfile(os.path.join(subdir_path, f)) and f.lower().endswith('.pdb')]

        if not pdb_files:
            print("  No PDB files found in this subdirectory. Skipping.")
            continue

        print(f"  Found {len(pdb_files)} PDB files. Extracting sequences...")

        alignment_data = [] # List to store (identifier, sequence) tuples
        processed_count = 0
        error_count = 0

        for pdb_filename in pdb_files:
            pdb_filepath = os.path.join(subdir_path, pdb_filename)
            identifier, sequence = extract_aligned_sequence(pdb_filepath, TARGET_SEQUENCE_LENGTH)

            if identifier and sequence:
                if len(sequence) == TARGET_SEQUENCE_LENGTH:
                    alignment_data.append((identifier, sequence))
                    processed_count += 1
                else:
                     # This check should ideally never fail if extract_aligned_sequence is correct
                     print(f"  Internal Error: Sequence length mismatch for {identifier} ({len(sequence)} != {TARGET_SEQUENCE_LENGTH}). Skipping.")
                     error_count += 1
            else:
                # Increment error count if extraction failed (warnings already printed)
                error_count += 1

        print(f"  Successfully extracted sequences for {processed_count} files.")
        if error_count > 0:
             print(f"  Failed to process or extract sequence from {error_count} files in this subdirectory.")

        # Write the A3M file if we have data
        if alignment_data:
            print(f"  Writing A3M alignment to: {a3m_output_file}")
            try:
                with open(a3m_output_file, 'w') as f_out:
                    # Write A3M header line (optional but common)
                    # f_out.write("# A3M\n")
                    for identifier, sequence in alignment_data:
                        f_out.write(f">{identifier}\n")
                        f_out.write(f"{sequence}\n")
                print(f"  Successfully wrote {len(alignment_data)} sequences.")
            except IOError as e:
                print(f"  Error writing A3M file {a3m_output_file}: {e}")
        elif processed_count == 0 and error_count > 0:
             print("  No valid sequences extracted due to errors, skipping A3M file generation.")
        else:
             # Case where subdir was empty or contained only non-PDB files
             print("  No valid sequence data extracted, skipping A3M file generation.")


    print("\n--- Script Finished ---")


# --- Command Line Argument Parsing ---
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Create A3M alignments from PDB segment files in categorized subdirectories.")
    parser.add_argument(
        "--inputdir",
        default=DEFAULT_INPUT_BASE_DIR,
        help=f"Base directory containing the PDB segment subdirectories (default: {DEFAULT_INPUT_BASE_DIR})"
    )
    parser.add_argument(
        "--subdirs",
        nargs='+', # Allows specifying multiple subdirectories
        default=DEFAULT_SUBDIRS,
        help=f"List of subdirectories to process (default: {' '.join(DEFAULT_SUBDIRS)})"
    )

    args = parser.parse_args()

    create_a3m_alignments(args.inputdir, args.subdirs)
