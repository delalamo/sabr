#!/usr/bin/env python3

import argparse
import os
import re

import pandas as pd

# Biopython is required: pip install biopython
try:
    import warnings

    from Bio.PDB import PDBIO, PDBParser, Select
    from Bio.PDB.PDBExceptions import PDBConstructionWarning
    from Bio.PDB.Polypeptide import is_aa

    # Suppress specific Biopython warnings if desired
    warnings.simplefilter('ignore', PDBConstructionWarning)
except ImportError:
    print("Error: Biopython library not found.")
    print("Please install it using: pip install biopython")
    exit(1)

# --- Configuration ---
# Default paths (can be overridden by command-line arguments)
DEFAULT_TSV_FILE_PATH = '/home/delalamo/sabdab_structures/sabdab_summary_all.tsv'
DEFAULT_PDB_DIR = '/home/delalamo/sabdab_structures/all_structures/imgt/'
DEFAULT_OUTPUT_BASE_DIR = './pdb_segments' # Output directories will be created here
RESIDUE_LIMIT = 128

# --- Helper Functions ---

def check_species(species_str, keywords):
    """
    Checks if a species string matches any of the keywords (case-insensitive).
    Handles potential non-string inputs gracefully.
    """
    if not isinstance(species_str, str):
        return False
    species_lower = species_str.lower()
    # Use word boundaries for more specific matching if needed, e.g., r'\bhuman\b'
    # For broader matching as requested (human, homo sapiens), simple 'in' is okay.
    return any(keyword.lower() in species_lower for keyword in keywords)

def is_valid_chain_id(chain_id):
    """
    Checks if a chain ID is valid (not NA, NaN, empty, etc.).
    """
    if pd.isna(chain_id):
        return False
    # Ensure it's treated as a string for comparison
    chain_id_str = str(chain_id).strip()
    if not chain_id_str or chain_id_str.upper() == 'NA':
        return False
    return True

class ResidueSelector(Select):
    """
    Biopython Select class to filter specific residues from a specific chain.
    Only accepts standard amino acid residues.
    """
    def __init__(self, chain_id, target_residues):
        self.chain_id = chain_id
        # Store target residue identifiers (tuple: (hetfield, resseq, icode))
        self.target_residue_ids = {res.id for res in target_residues}

    def accept_chain(self, chain):
        # Process only the target chain
        return chain.id == self.chain_id

    def accept_residue(self, residue):
        # Accept residue only if it's in our target list AND is a standard amino acid
        # residue.id[0] == ' ' checks for hetfield, ' ' means standard residue
        return residue.id in self.target_residue_ids and residue.id[0] == ' ' and is_aa(residue, standard=True)

    def accept_atom(self, atom):
        # Accept all atoms of the accepted residues
        return True

# --- Main Processing Logic ---

def process_pdb_files(tsv_file_path, pdb_dir, output_base_dir):
    """
    Reads the TSV, identifies target chains, extracts residues from PDBs,
    and saves them to categorized directories.
    """
    print("--- PDB Chain Segment Extractor ---")
    print(f"TSV file: {tsv_file_path}")
    print(f"PDB directory: {pdb_dir}")
    print(f"Output base directory: {output_base_dir}")
    print(f"Residue limit per chain: {RESIDUE_LIMIT}")

    # Define output directories relative to the output base directory
    output_dirs = {
        'human_heavy': os.path.join(output_base_dir, 'human_heavy'),
        'human_light': os.path.join(output_base_dir, 'human_light'),
        'mouse_heavy': os.path.join(output_base_dir, 'mouse_heavy'),
        'mouse_light': os.path.join(output_base_dir, 'mouse_light'),
        'camelid_heavy': os.path.join(output_base_dir, 'camelid_heavy')
    }

    # Create output directories if they don't exist
    print("\nCreating output directories...")
    for dir_path in output_dirs.values():
        try:
            os.makedirs(dir_path, exist_ok=True)
            print(f"  Ensured directory exists: {dir_path}")
        except OSError as e:
            print(f"Error creating directory {dir_path}: {e}")
            return # Stop if we can't create output directories

    # --- Read TSV File ---
    print(f"\nReading TSV file: {tsv_file_path}...")
    try:
        # Attempt to read with tab separator first
        try:
            df = pd.read_csv(tsv_file_path, sep='\t', on_bad_lines='warn')
        except pd.errors.ParserError:
            print("  Tab separator failed, trying whitespace separator...")
            # If tab fails, try with flexible whitespace separator
            df = pd.read_csv(tsv_file_path, sep='\s+', engine='python', on_bad_lines='warn')

        # Clean up column names (remove leading/trailing spaces)
        df.columns = [col.strip() for col in df.columns]
        print(f"Successfully read {len(df)} rows.")

        # Verify required columns exist
        required_cols = ['pdb', 'Hchain', 'Lchain', 'heavy_species', 'light_species']
        missing_cols = [col for col in required_cols if col not in df.columns]
        if missing_cols:
            print(f"\nError: Missing required columns in TSV: {missing_cols}")
            print(f"Found columns: {df.columns.tolist()}")
            return
        print("Required columns found.")

    except FileNotFoundError:
        print(f"\nError: TSV file not found at {tsv_file_path}")
        return
    except Exception as e:
        print(f"\nError reading or parsing TSV file: {e}")
        return

    # --- Initialize Counters and PDB Parser ---
    processed_rows = 0
    target_chain_tasks = [] # List to store (pdb_code, chain_id, target_dir_key) tuples
    parser = PDBParser(QUIET=True) # QUIET=True suppresses many warnings

    # --- Identify Target Chains from TSV ---
    print("\nIdentifying target chains from TSV...")
    for index, row in df.iterrows():
        processed_rows += 1
        pdb_code = str(row['pdb']).lower().strip()
        if not pdb_code or len(pdb_code) != 4:
            # print(f"  Skipping row {index+2}: Invalid or missing PDB code '{row['pdb']}'")
            continue

        # Check Heavy Chain
        h_chain_id_raw = row['Hchain']
        h_species = row['heavy_species']
        if is_valid_chain_id(h_chain_id_raw):
            h_chain_id = str(h_chain_id_raw).strip() # Ensure it's a clean string
            target_dir_key_h = None
            if check_species(h_species, ['homo sapiens', 'human']):
                target_dir_key_h = 'human_heavy'
            elif check_species(h_species, ['mus musculus', 'mouse']):
                target_dir_key_h = 'mouse_heavy'
            # Broader check for camelids
            elif check_species(h_species, ['camel', 'lama', 'vicugna', 'alpaca', 'dromaderius']):
                target_dir_key_h = 'camelid_heavy'

            if target_dir_key_h:
                target_chain_tasks.append((pdb_code, h_chain_id, target_dir_key_h))

        # Check Light Chain
        l_chain_id_raw = row['Lchain']
        l_species = row['light_species']
        if is_valid_chain_id(l_chain_id_raw):
            l_chain_id = str(l_chain_id_raw).strip() # Ensure it's a clean string
            target_dir_key_l = None
            if check_species(l_species, ['homo sapiens', 'human']):
                target_dir_key_l = 'human_light'
            elif check_species(l_species, ['mus musculus', 'mouse']):
                target_dir_key_l = 'mouse_light'

            if target_dir_key_l:
                target_chain_tasks.append((pdb_code, l_chain_id, target_dir_key_l))

    print(f"Identified {len(target_chain_tasks)} potential chain extraction tasks from {processed_rows} processed rows.")

    # --- Process PDB Files for Identified Chains ---
    print("\nProcessing PDB files and extracting segments...")
    io = PDBIO()
    saved_files_count = {key: 0 for key in output_dirs}
    skipped_missing_pdb = 0
    skipped_missing_chain = 0
    skipped_no_aa = 0
    processed_short_chain = 0 # Chains processed but shorter than limit
    error_parsing_pdb = 0
    processed_tasks = 0

    # Use a set to avoid processing the same PDB file multiple times if not needed,
    # but here we need to process per chain task.
    unique_tasks = set(target_chain_tasks) # Remove duplicate tasks if any
    print(f"Processing {len(unique_tasks)} unique (PDB code, chain ID, category) tasks.")

    for pdb_code, chain_id, target_dir_key in unique_tasks:
        processed_tasks += 1
        if processed_tasks % 100 == 0: # Print progress update
             print(f"  Processed {processed_tasks}/{len(unique_tasks)} tasks...")

        pdb_file = os.path.join(pdb_dir, f"{pdb_code}.pdb")

        if not os.path.exists(pdb_file):
            # print(f"  Skipping task ({pdb_code}, {chain_id}): PDB file not found at {pdb_file}")
            skipped_missing_pdb += 1
            continue

        try:
            # Parse the structure
            structure = parser.get_structure(pdb_code, pdb_file)

            # Access the first model (assuming single model PDB, common case)
            if len(structure) == 0:
                 # print(f"  Skipping task ({pdb_code}, {chain_id}): No models found in PDB.")
                 error_parsing_pdb += 1 # Count as parsing error
                 continue
            model = structure[0]

            # Check if the target chain exists in the model
            if chain_id not in model:
                # print(f"  Skipping task ({pdb_code}, {chain_id}): Chain '{chain_id}' not found in model 0.")
                skipped_missing_chain += 1
                continue

            chain = model[chain_id]
            residues_to_save = []
            aa_count = 0

            # Iterate through residues and collect the first RESIDUE_LIMIT standard amino acids
            for residue in chain:
                # Check if it's a standard amino acid (hetfield==' ' and is_aa)
                if residue.id[0] == ' ' and is_aa(residue, standard=True):
                    residues_to_save.append(residue)
                    aa_count += 1
                    if aa_count >= RESIDUE_LIMIT:
                        break # Stop once we reach the limit

            # Proceed if we found any standard amino acid residues
            if aa_count > 0:
                # If chain has fewer standard AAs than the limit, note it
                if aa_count < RESIDUE_LIMIT:
                    # print(f"  Note: Chain {chain_id} in {pdb_code} has only {aa_count} standard AA residues (limit {RESIDUE_LIMIT}). Saving truncated segment.")
                    processed_short_chain += 1

                # Define the output filename
                output_filename = os.path.join(output_dirs[target_dir_key], f"{pdb_code}_{chain_id}_1-{aa_count}.pdb")

                # Save the selected residues using the custom selector
                io.set_structure(structure) # Provide the whole structure context
                io.save(output_filename, ResidueSelector(chain_id, residues_to_save))
                saved_files_count[target_dir_key] += 1

            else:
                # print(f"  Skipping task ({pdb_code}, {chain_id}): No standard amino acid residues found in chain.")
                skipped_no_aa += 1

        except Exception as e:
            # Catch potential errors during PDB parsing or processing
            # print(f"  Error processing PDB task ({pdb_code}, {chain_id}): {e}")
            error_parsing_pdb += 1

    # --- Final Summary ---
    print("\n--- Processing Summary ---")
    print(f"Total rows processed from TSV: {processed_rows}")
    print(f"Total unique extraction tasks identified: {len(unique_tasks)}")
    print(f"Tasks skipped due to missing PDB file: {skipped_missing_pdb}")
    print(f"Tasks skipped due to missing chain in PDB model: {skipped_missing_chain}")
    print(f"Tasks skipped due to no standard AAs found in chain: {skipped_no_aa}")
    print(f"Tasks skipped due to PDB parsing/processing errors: {error_parsing_pdb}")
    print(f"Chains processed but shorter than {RESIDUE_LIMIT} residues: {processed_short_chain}")

    print("\nFiles saved:")
    total_saved = 0
    for dir_key, count in saved_files_count.items():
        print(f"  - {output_dirs[dir_key]}: {count} files")
        total_saved += count

    print(f"\nTotal PDB segment files saved: {total_saved}")
    print("--- Script Finished ---")


# --- Command Line Argument Parsing ---
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Extract N-terminal segments of specific antibody chains from PDB files based on a TSV summary.")
    parser.add_argument(
        "--tsv",
        default=DEFAULT_TSV_FILE_PATH,
        help=f"Path to the input TSV summary file (default: {DEFAULT_TSV_FILE_PATH})"
    )
    parser.add_argument(
        "--pdbdir",
        default=DEFAULT_PDB_DIR,
        help=f"Path to the directory containing PDB files (default: {DEFAULT_PDB_DIR})"
    )
    parser.add_argument(
        "--outdir",
        default=DEFAULT_OUTPUT_BASE_DIR,
        help=f"Base directory where output category folders will be created (default: {DEFAULT_OUTPUT_BASE_DIR})"
    )

    args = parser.parse_args()

    process_pdb_files(args.tsv, args.pdbdir, args.outdir)
