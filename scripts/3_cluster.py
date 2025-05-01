#!/usr/bin/env python3

import argparse
import math
import os
from collections import defaultdict

# Biopython is required for AlignIO: pip install biopython
try:
    from Bio import AlignIO
    from Bio.Seq import Seq
    # Removed: from Bio.SubsMat import MatrixInfo
    from Bio.SeqRecord import SeqRecord
except ImportError:
    print("Error: Biopython library not found (required for AlignIO).")
    print("Please install it using: pip install biopython")
    exit(1)
except Exception as e:
    print(f"Error importing Biopython components: {e}")
    exit(1)

# --- Configuration ---
# Default input/output directory
DEFAULT_BASE_DIR = './pdb_segments'
# Default subdirectories (used to find .a3m files)
DEFAULT_SUBDIRS = [
    'human_heavy', 'human_light', 'mouse_heavy',
    'mouse_light', 'camelid_heavy'
]
# Default clustering threshold for normalized BLOSUM62 score
DEFAULT_THRESHOLD = 0.7
# Gap penalty for residue-gap alignment in score calculation
GAP_PENALTY = -4
# Character for unknown residues
UNKNOWN_CHAR = 'X'
# Gap character
GAP_CHAR = '-'

# --- Define BLOSUM62 Matrix Directly (User Provided) ---
# This dictionary contains the BLOSUM62 scores.
# The lookup function below handles key order (e.g., ('A', 'R') vs ('R', 'A')).
BLOSUM62_MATRIX = {
    ('W', 'F'): 1, ('L', 'R'): -2, ('S', 'P'): -1, ('V', 'T'): 0,
    ('Q', 'Q'): 5, ('N', 'A'): -2, ('Z', 'Y'): -2, ('W', 'R'): -3,
    ('Q', 'A'): -1, ('S', 'D'): 0, ('H', 'H'): 8, ('S', 'H'): -1,
    ('H', 'D'): -1, ('L', 'N'): -3, ('W', 'A'): -3, ('Y', 'M'): -1,
    ('G', 'R'): -2, ('Y', 'I'): -1, ('Y', 'E'): -2, ('B', 'Y'): -3,
    ('Y', 'A'): -2, ('V', 'D'): -3, ('B', 'S'): 0, ('Y', 'Y'): 7,
    ('G', 'N'): 0, ('E', 'C'): -4, ('Y', 'Q'): -1, ('Z', 'Z'): 4,
    ('V', 'A'): 0, ('C', 'C'): 9, ('M', 'R'): -1, ('V', 'E'): -2,
    ('T', 'N'): 0, ('P', 'P'): 7, ('V', 'I'): 3, ('V', 'S'): -2,
    ('Z', 'P'): -1, ('V', 'M'): 1, ('T', 'F'): -2, ('V', 'Q'): -2,
    ('K', 'K'): 5, ('P', 'D'): -1, ('I', 'H'): -3, ('I', 'D'): -3,
    ('T', 'R'): -1, ('P', 'L'): -3, ('K', 'G'): -2, ('M', 'N'): -2,
    ('P', 'H'): -2, ('F', 'Q'): -3, ('Z', 'G'): -2, ('X', 'L'): -1,
    ('T', 'M'): -1, ('Z', 'C'): -3, ('X', 'H'): -1, ('D', 'R'): -2,
    ('B', 'W'): -4, ('X', 'D'): -1, ('Z', 'K'): 1, ('F', 'A'): -2,
    ('Z', 'W'): -3, ('F', 'E'): -3, ('D', 'N'): 1, ('B', 'K'): 0,
    ('X', 'X'): -1, ('F', 'I'): 0, ('B', 'G'): -1, ('X', 'T'): 0,
    ('F', 'M'): 0, ('B', 'C'): -3, ('Z', 'I'): -3, ('Z', 'V'): -2,
    ('S', 'S'): 4, ('L', 'Q'): -2, ('W', 'E'): -3, ('Q', 'R'): 1,
    ('N', 'N'): 6, ('W', 'M'): -1, ('Q', 'C'): -3, ('W', 'I'): -3,
    ('S', 'C'): -1, ('L', 'A'): -1, ('S', 'G'): 0, ('L', 'E'): -3,
    ('W', 'Q'): -2, ('H', 'G'): -2, ('S', 'K'): 0, ('Q', 'N'): 0,
    ('N', 'R'): 0, ('H', 'C'): -3, ('Y', 'N'): -2, ('G', 'Q'): -2,
    ('Y', 'F'): 3, ('C', 'A'): 0, ('V', 'L'): 1, ('G', 'E'): -2,
    ('G', 'A'): 0, ('K', 'R'): 2, ('E', 'D'): 2, ('Y', 'R'): -2,
    ('M', 'Q'): 0, ('T', 'I'): -1, ('C', 'D'): -3, ('V', 'F'): -1,
    ('T', 'A'): 0, ('T', 'P'): -1, ('B', 'P'): -2, ('T', 'E'): -1,
    ('V', 'N'): -3, ('P', 'G'): -2, ('M', 'A'): -1, ('K', 'H'): -1,
    ('V', 'R'): -3, ('P', 'C'): -3, ('M', 'E'): -2, ('K', 'L'): -2,
    ('V', 'V'): 4, ('M', 'I'): 1, ('T', 'Q'): -1, ('I', 'G'): -4,
    ('P', 'K'): -1, ('M', 'M'): 5, ('K', 'D'): -1, ('I', 'C'): -1,
    ('Z', 'D'): 1, ('F', 'R'): -3, ('X', 'K'): -1, ('Q', 'D'): 0,
    ('X', 'G'): -1, ('Z', 'L'): -3, ('X', 'C'): -2, ('Z', 'H'): 0,
    ('B', 'L'): -4, ('B', 'H'): 0, ('F', 'F'): 6, ('X', 'W'): -2,
    ('B', 'D'): 4, ('D', 'A'): -2, ('S', 'L'): -2, ('X', 'S'): 0,
    ('F', 'N'): -3, ('S', 'R'): -1, ('W', 'D'): -4, ('V', 'Y'): -1,
    ('W', 'L'): -2, ('H', 'R'): 0, ('W', 'H'): -2, ('H', 'N'): 1,
    ('W', 'T'): -2, ('T', 'T'): 5, ('S', 'F'): -2, ('W', 'P'): -4,
    ('L', 'D'): -4, ('B', 'I'): -3, ('L', 'H'): -3, ('S', 'N'): 1,
    ('B', 'T'): -1, ('L', 'L'): 4, ('Y', 'K'): -2, ('E', 'Q'): 2,
    ('Y', 'G'): -3, ('Z', 'S'): 0, ('Y', 'C'): -2, ('G', 'D'): -1,
    ('B', 'V'): -3, ('E', 'A'): -1, ('Y', 'W'): 2, ('E', 'E'): 5,
    ('Y', 'S'): -2, ('C', 'N'): -3, ('V', 'C'): -1, ('T', 'H'): -2,
    ('P', 'R'): -2, ('V', 'G'): -3, ('T', 'L'): -1, ('V', 'K'): -2,
    ('K', 'Q'): 1, ('R', 'A'): -1, ('I', 'R'): -3, ('T', 'D'): -1,
    ('P', 'F'): -4, ('I', 'N'): -3, ('K', 'I'): -3, ('M', 'D'): -3,
    ('V', 'W'): -3, ('W', 'W'): 11, ('M', 'H'): -2, ('P', 'N'): -2,
    ('K', 'A'): -1, ('M', 'L'): 2, ('K', 'E'): 1, ('Z', 'E'): 4,
    ('X', 'N'): -1, ('Z', 'A'): -1, ('Z', 'M'): -1, ('X', 'F'): -1,
    ('K', 'C'): -3, ('B', 'Q'): 0, ('X', 'B'): -1, ('B', 'M'): -3,
    ('F', 'C'): -2, ('Z', 'Q'): 3, ('X', 'Z'): -1, ('F', 'G'): -3,
    ('B', 'E'): 1, ('X', 'V'): -1, ('F', 'K'): -3, ('B', 'A'): -2,
    ('X', 'R'): -1, ('D', 'D'): 6, ('W', 'G'): -2, ('Z', 'F'): -3,
    ('S', 'Q'): 0, ('W', 'C'): -2, ('W', 'K'): -3, ('H', 'Q'): 0,
    ('L', 'C'): -1, ('W', 'N'): -4, ('S', 'A'): 1, ('L', 'G'): -4,
    ('W', 'S'): -3, ('S', 'E'): 0, ('H', 'E'): 0, ('S', 'I'): -2,
    ('H', 'A'): -2, ('S', 'M'): -1, ('Y', 'L'): -1, ('Y', 'H'): 2,
    ('Y', 'D'): -3, ('E', 'R'): 0, ('X', 'P'): -2, ('G', 'G'): 6,
    ('G', 'C'): -3, ('E', 'N'): 0, ('Y', 'T'): -2, ('Y', 'P'): -3,
    ('T', 'K'): -1, ('A', 'A'): 4, ('P', 'Q'): -1, ('T', 'C'): -1,
    ('V', 'H'): -3, ('T', 'G'): -2, ('I', 'Q'): -3, ('Z', 'T'): -1,
    ('C', 'R'): -3, ('V', 'P'): -2, ('P', 'E'): -1, ('M', 'C'): -1,
    ('K', 'N'): 0, ('I', 'I'): 4, ('P', 'A'): -1, ('M', 'G'): -3,
    ('T', 'S'): 1, ('I', 'E'): -3, ('P', 'M'): -2, ('M', 'K'): -1,
    ('I', 'A'): -1, ('P', 'I'): -3, ('R', 'R'): 5, ('X', 'M'): -1,
    ('L', 'I'): 2, ('X', 'I'): -1, ('Z', 'B'): 1, ('X', 'E'): -1,
    ('Z', 'N'): 0, ('X', 'A'): 0, ('B', 'R'): -1, ('B', 'N'): 3,
    ('F', 'D'): -3, ('X', 'Y'): -1, ('Z', 'R'): 0, ('F', 'H'): -1,
    ('B', 'F'): -3, ('F', 'L'): 0, ('X', 'Q'): -1, ('B', 'B'): 4
}

# --- Helper Functions ---

def get_blosum_score(res1, res2, matrix):
    """
    Looks up the BLOSUM62 score for a pair of residues using the provided matrix.
    Handles key order (e.g., (A, R) vs (R, A)) and missing pairs.
    """
    # Ensure consistent key order by sorting the residue pair
    pair = tuple(sorted((res1, res2)))
    # Return the score from the matrix, defaulting to 0 if the pair is not found
    # (e.g., for non-standard residues if they weren't filtered earlier)
    return matrix.get(pair, 0)

def calculate_blosum_score(seq1, seq2, matrix, gap_penalty=GAP_PENALTY):
    """
    Calculates the raw BLOSUM62 score between two aligned sequences.
    Handles gaps and unknown ('X') residues using the provided matrix.
    """
    score = 0
    if len(seq1) != len(seq2):
        raise ValueError("Sequences must have the same length for aligned score calculation.")

    alignment_length = len(seq1)
    for i in range(alignment_length):
        res1 = seq1[i]
        res2 = seq2[i]
        pair_score = 0 # Default score

        if res1 == GAP_CHAR and res2 == GAP_CHAR:
            pair_score = 0  # Gap-Gap score is 0
        elif res1 == GAP_CHAR or res2 == GAP_CHAR:
            pair_score = gap_penalty  # Residue-Gap incurs penalty
        elif res1 == UNKNOWN_CHAR or res2 == UNKNOWN_CHAR:
            pair_score = 0 # Treat 'X' neutrally against anything
        else:
            # Standard residue pair - look up in the provided BLOSUM matrix
            pair_score = get_blosum_score(res1, res2, matrix)

        score += pair_score
    return score

def calculate_normalized_score(id1, id2, sequences, matrix, self_scores, gap_penalty=GAP_PENALTY):
    """
    Calculates the normalized BLOSUM62 score between two sequences,
    using pre-calculated self-scores for normalization.
    Normalization is score / min(self_score1, self_score2).
    Uses the provided matrix.
    """
    seq1 = sequences[id1]
    seq2 = sequences[id2]

    # Retrieve precomputed self-scores
    self_score1 = self_scores[id1]
    self_score2 = self_scores[id2]

    # Determine the denominator for normalization
    denominator = min(self_score1, self_score2)

    # Handle cases where self-similarity isn't positive (prevents division errors)
    # This can happen if a sequence is all gaps or has very low self-score.
    if denominator <= 0:
        # If sequences cannot even score positively against themselves,
        # similarity to others is considered zero.
        # Check if pairwise score is also non-positive to avoid 0/0 -> 1.0 issues
        pairwise_score_check = calculate_blosum_score(seq1, seq2, matrix, gap_penalty)
        return 0.0 if pairwise_score_check <= 0 else 1.0 # Or handle as error/special case

    # Calculate the pairwise score between the two different sequences
    pairwise_score = calculate_blosum_score(seq1, seq2, matrix, gap_penalty)

    # Normalize the score
    normalized = pairwise_score / denominator

    # Clamp the score to be at least 0.0
    return max(0.0, normalized)


def greedy_cluster(sequences, matrix, threshold, gap_penalty=GAP_PENALTY):
    """
    Performs greedy clustering based on normalized BLOSUM62 score.
    Uses the provided matrix.

    Args:
        sequences (dict): Dictionary {identifier: sequence_string}.
        matrix (dict): The BLOSUM62 substitution matrix dictionary.
        threshold (float): The normalized score threshold for clustering.
        gap_penalty (int): Gap penalty used in score calculation.

    Returns:
        list: A list of clusters, where each cluster is a list of sequence identifiers.
    """
    if not sequences:
        return []

    print(f"  Starting greedy clustering with threshold {threshold:.2f}...")
    num_sequences = len(sequences)
    # Sort IDs for deterministic clustering (optional but good practice)
    seq_ids = sorted(list(sequences.keys()))

    # --- Precompute self-scores ---
    print("    Precomputing self-scores...")
    self_scores = {}
    for seq_id in seq_ids:
        seq = sequences[seq_id]
        # Calculate self-score using the provided matrix
        self_scores[seq_id] = calculate_blosum_score(seq, seq, matrix, gap_penalty)
        # print(f"      Self-score for {seq_id}: {self_scores[seq_id]}") # Debugging

    clustered_ids = set() # Keep track of IDs already assigned to a cluster
    clusters = [] # List to store the final clusters

    print("    Running clustering iterations...")
    processed_count = 0
    # Iterate through sequences to potentially form new cluster centroids
    for i, centroid_id in enumerate(seq_ids):
        if centroid_id in clustered_ids:
            continue # Skip if already clustered

        processed_count += 1
        if processed_count % 100 == 0 or processed_count == num_sequences:
             print(f"      Processed {processed_count}/{num_sequences} potential centroids...")

        # Start a new cluster with the current sequence as the centroid
        current_cluster = [centroid_id]
        clustered_ids.add(centroid_id)

        # Compare this centroid to all other *unclustered* sequences
        for j in range(i + 1, len(seq_ids)): # Only need to check subsequent IDs
            candidate_id = seq_ids[j]
            if candidate_id in clustered_ids:
                continue # Skip if already clustered

            # Calculate normalized similarity score using the provided matrix
            norm_score = calculate_normalized_score(
                centroid_id, candidate_id, sequences, matrix, self_scores, gap_penalty
            )

            # If score meets threshold, add candidate to the current cluster
            if norm_score >= threshold:
                current_cluster.append(candidate_id)
                clustered_ids.add(candidate_id)

        clusters.append(current_cluster) # Add the completed cluster

    print(f"  Clustering finished. Found {len(clusters)} clusters.")
    return clusters

# --- Main Processing Logic ---

def process_alignments(base_dir, subdirs_to_process, threshold):
    """
    Reads A3M files, performs clustering using the internal BLOSUM62 matrix,
    and writes cluster files.
    """
    print("--- Alignment-Based BLOSUM62 Clustering Script ---")
    print(f"Base directory: {base_dir}")
    print(f"Processing subdirectories: {', '.join(subdirs_to_process)}")
    print(f"Normalized BLOSUM62 score threshold: {threshold}")
    print(f"Gap Penalty: {GAP_PENALTY}")

    # Use the internally defined BLOSUM62 matrix
    matrix = BLOSUM62_MATRIX

    if not os.path.isdir(base_dir):
        print(f"Error: Base directory not found: {base_dir}")
        return

    # Process each specified subdirectory type
    for subdir_name in subdirs_to_process:
        a3m_input_file = os.path.join(base_dir, f"{subdir_name}.a3m")
        cluster_output_file = os.path.join(base_dir, f"{subdir_name}.clusters")

        print(f"\nProcessing alignment file: {a3m_input_file}...")

        if not os.path.isfile(a3m_input_file):
            print(f"  Warning: Alignment file '{a3m_input_file}' not found. Skipping.")
            continue

        # --- Read Alignment File ---
        sequences = {}
        try:
            # Use AlignIO to read the aligned FASTA/A3M format
            alignment = AlignIO.read(a3m_input_file, "fasta")
            for record in alignment:
                # Ensure sequence is uppercase and store it
                sequences[record.id] = str(record.seq).upper()
            print(f"  Read {len(sequences)} sequences from alignment.")
            if not sequences:
                 print("  Alignment file is empty. Skipping clustering.")
                 continue

        except ValueError as e:
             # Catch specific AlignIO errors like empty file or wrong format
             print(f"  Error reading alignment file {a3m_input_file}: {e}. Is it a valid FASTA/A3M alignment?")
             continue
        except Exception as e:
            print(f"  Unexpected error reading alignment file {a3m_input_file}: {e}")
            continue # Skip to the next file

        # --- Perform Clustering ---
        # Pass the internal matrix dictionary to the clustering function
        clusters = greedy_cluster(sequences, matrix, threshold, GAP_PENALTY)

        # --- Write Cluster Output File ---
        if clusters:
            print(f"  Writing cluster results to: {cluster_output_file}")
            try:
                with open(cluster_output_file, 'w') as f_out:
                    for i, cluster in enumerate(clusters):
                        # Write space-separated identifiers for each cluster
                        f_out.write(" ".join(cluster) + "\n")
                print(f"  Successfully wrote {len(clusters)} clusters.")
            except IOError as e:
                print(f"  Error writing cluster file {cluster_output_file}: {e}")
        else:
            print("  No clusters generated (or input was empty/unreadable). No output file written.")

    print("\n--- Script Finished ---")


# --- Command Line Argument Parsing ---
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Cluster sequences in A3M alignments based on normalized BLOSUM62 score using greedy clustering.")
    parser.add_argument(
        "--inputdir",
        default=DEFAULT_BASE_DIR,
        help=f"Base directory containing the .a3m alignment files and where .clusters files will be written (default: {DEFAULT_BASE_DIR})"
    )
    parser.add_argument(
        "--subdirs",
        nargs='+', # Allows specifying multiple subdirectory names (without path)
        default=DEFAULT_SUBDIRS,
        help=f"List of subdirectory names used to identify .a3m files (e.g., human_heavy human_light) (default: {' '.join(DEFAULT_SUBDIRS)})"
    )
    parser.add_argument(
        "--threshold",
        type=float,
        default=DEFAULT_THRESHOLD,
        help=f"Normalized BLOSUM62 score threshold for clustering (default: {DEFAULT_THRESHOLD})"
    )

    args = parser.parse_args()

    # Validate threshold range slightly differently, as normalized scores might exceed 1
    if args.threshold < 0:
         print(f"Error: Threshold cannot be negative. Received: {args.threshold}")
         exit(1)
    elif args.threshold > 1.5: # Arbitrary upper bound check to catch likely typos
         print(f"Warning: Threshold {args.threshold} seems high for normalized score. Ensure this is intended.")


    process_alignments(args.inputdir, args.subdirs, args.threshold)
