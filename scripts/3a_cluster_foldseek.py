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
DEFAULT_THRESHOLD = 0.8
# Gap penalty for residue-gap alignment in score calculation
GAP_PENALTY = -1
GAP_OPEN_PENALTY = -10
# Character for unknown residues
UNKNOWN_CHAR = 'X'
# Gap character
GAP_CHAR = '-'

# --- Define Foldseek Matrix Directly (User Provided) ---
# This dictionary contains the BLOSUM62 scores.
# The lookup function below handles key order (e.g., ('A', 'R') vs ('R', 'A')).


BLOSUM62_MATRIX = {
        ("A", "A"): 6,  ("A", "C"): -3,  ("A", "D"):  1, ("A", "E"):  2, ("A", "F"):  3, ("A", "G"): -2, ("A", "H"): -2, ("A", "I"): -7, ("A", "K"): -3, ("A", "L"): -3, ("A", "M"):-10, ("A", "N"): -5, ("A", "P"): -1, ("A", "Q"):  1, ("A", "R"): -4, ("A", "S"): -7, ("A", "T"): -5, ("A", "V"): -6, ("A", "W"):  0, ("A", "Y"): -2, ("A", "X"): 0,
        ("C", "A"): -3, ("C", "C"):  6,  ("C", "D"): -2, ("C", "E"): -8, ("C", "F"): -5, ("C", "G"): -4, ("C", "H"): -4, ("C", "I"):-12, ("C", "K"):-13, ("C", "L"):  1, ("C", "M"):-14, ("C", "N"):  0, ("C", "P"):  0, ("C", "Q"):  1, ("C", "R"): -1, ("C", "S"):  0, ("C", "T"): -8, ("C", "V"):  1, ("C", "W"): -7, ("C", "Y"): -9, ("C", "X"): 0,
        ("D", "A"):  1, ("D", "C"): -2,  ("D", "D"):  4, ("D", "E"): -3, ("D", "F"):  0, ("D", "G"):  1, ("D", "H"):  1, ("D", "I"): -3, ("D", "K"): -5, ("D", "L"): -4, ("D", "M"): -5, ("D", "N"): -2, ("D", "P"):  1, ("D", "Q"): -1, ("D", "R"): -1, ("D", "S"): -4, ("D", "T"): -2, ("D", "V"): -3, ("D", "W"): -2, ("D", "Y"): -2, ("D", "X"): 0,
        ("E", "A"):  2, ("E", "C"): -8,  ("E", "D"): -3, ("E", "E"):  9, ("E", "F"): -2, ("E", "G"): -7, ("E", "H"): -4, ("E", "I"):-12, ("E", "K"):-10, ("E", "L"): -7, ("E", "M"):-17, ("E", "N"): -8, ("E", "P"): -6, ("E", "Q"): -3, ("E", "R"): -8, ("E", "S"):-10, ("E", "T"):-10, ("E", "V"):-13, ("E", "W"): -6, ("E", "Y"): -3, ("E", "X"): 0,
        ("F", "A"):  3, ("F", "C"): -5,  ("F", "D"):  0, ("F", "E"): -2, ("F", "F"):  7, ("F", "G"): -3, ("F", "H"): -3, ("F", "I"): -5, ("F", "K"):  1, ("F", "L"): -3, ("F", "M"): -9, ("F", "N"): -5, ("F", "P"): -2, ("F", "Q"):  2, ("F", "R"): -5, ("F", "S"): -8, ("F", "T"): -3, ("F", "V"): -7, ("F", "W"):  4, ("F", "Y"): -4, ("F", "X"): 0,
        ("G", "A"): -2, ("G", "C"): -4,  ("G", "D"):  1, ("G", "E"): -7, ("G", "F"): -3, ("G", "G"):  6, ("G", "H"):  3, ("G", "I"):  0, ("G", "K"): -7, ("G", "L"): -7, ("G", "M"): -1, ("G", "N"): -2, ("G", "P"): -2, ("G", "Q"): -4, ("G", "R"):  3, ("G", "S"): -3, ("G", "T"):  4, ("G", "V"): -6, ("G", "W"): -4, ("G", "Y"): -2, ("G", "X"): 0,
        ("H", "A"): -2, ("H", "C"): -4,  ("H", "D"):  1, ("H", "E"): -4, ("H", "F"): -3, ("H", "G"):  3, ("H", "H"):  6, ("H", "I"): -4, ("H", "K"): -7, ("H", "L"): -6, ("H", "M"): -6, ("H", "N"):  0, ("H", "P"): -1, ("H", "Q"): -3, ("H", "R"):  1, ("H", "S"): -3, ("H", "T"): -1, ("H", "V"): -5, ("H", "W"): -5, ("H", "Y"):  3, ("H", "X"): 0,
        ("I", "A"): -7, ("I", "C"):-12,  ("I", "D"): -3, ("I", "E"):-12, ("I", "F"): -5, ("I", "G"):  0, ("I", "H"): -4, ("I", "I"):  8, ("I", "K"): -5, ("I", "L"):-11, ("I", "M"):  7, ("I", "N"): -7, ("I", "P"): -6, ("I", "Q"): -6, ("I", "R"): -3, ("I", "S"): -9, ("I", "T"):  6, ("I", "V"):-12, ("I", "W"): -5, ("I", "Y"): -8, ("I", "X"): 0,
        ("K", "A"): -3, ("K", "C"):-13,  ("K", "D"): -5, ("K", "E"):-10, ("K", "F"):  1, ("K", "G"): -7, ("K", "H"): -7, ("K", "I"): -5, ("K", "K"):  9, ("K", "L"):-11, ("K", "M"): -8, ("K", "N"):-12, ("K", "P"): -6, ("K", "Q"): -5, ("K", "R"): -9, ("K", "S"):-14, ("K", "T"): -5, ("K", "V"):-15, ("K", "W"):  5, ("K", "Y"): -8, ("K", "X"): 0,
        ("L", "A"): -3, ("L", "C"):  1,  ("L", "D"): -4, ("L", "E"): -7, ("L", "F"): -3, ("L", "G"): -7, ("L", "H"): -6, ("L", "I"):-11, ("L", "K"):-11, ("L", "L"):  6, ("L", "M"):-16, ("L", "N"): -3, ("L", "P"): -2, ("L", "Q"):  2, ("L", "R"): -4, ("L", "S"): -4, ("L", "T"): -9, ("L", "V"):  0, ("L", "W"): -8, ("L", "Y"): -9, ("L", "X"): 0,
        ("M", "A"):-10, ("M", "C"):-14,  ("M", "D"): -5, ("M", "E"):-17, ("M", "F"): -9, ("M", "G"): -1, ("M", "H"): -6, ("M", "I"):  7, ("M", "K"): -8, ("M", "L"):-16, ("M", "M"): 10, ("M", "N"): -9, ("M", "P"): -9, ("M", "Q"):-10, ("M", "R"): -5, ("M", "S"):-10, ("M", "T"):  3, ("M", "V"):-16, ("M", "W"): -6, ("M", "Y"): -9, ("M", "X"): 0,
        ("N", "A"): -5, ("N", "C"):  0,  ("N", "D"): -2, ("N", "E"): -8, ("N", "F"): -5, ("N", "G"): -2, ("N", "H"):  0, ("N", "I"): -7, ("N", "K"):-12, ("N", "L"): -3, ("N", "M"): -9, ("N", "N"):  7, ("N", "P"):  0, ("N", "Q"): -2, ("N", "R"):  2, ("N", "S"):  3, ("N", "T"): -4, ("N", "V"):  0, ("N", "W"): -8, ("N", "Y"): -5, ("N", "X"): 0,
        ("P", "A"): -1, ("P", "C"):  0,  ("P", "D"):  1, ("P", "E"): -6, ("P", "F"): -2, ("P", "G"): -2, ("P", "H"): -1, ("P", "I"): -6, ("P", "K"): -6, ("P", "L"): -2, ("P", "M"): -9, ("P", "N"):  0, ("P", "P"):  4, ("P", "Q"):  0, ("P", "R"):  0, ("P", "S"): -2, ("P", "T"): -4, ("P", "V"):  0, ("P", "W"): -4, ("P", "Y"): -5, ("P", "X"): 0,
        ("Q", "A"):  1, ("Q", "C"):  1,  ("Q", "D"): -1, ("Q", "E"): -3, ("Q", "F"):  2, ("Q", "G"): -4, ("Q", "H"): -3, ("Q", "I"): -6, ("Q", "K"): -5, ("Q", "L"):  2, ("Q", "M"):-10, ("Q", "N"): -2, ("Q", "P"):  0, ("Q", "Q"):  5, ("Q", "R"): -2, ("Q", "S"): -4, ("Q", "T"): -5, ("Q", "V"): -1, ("Q", "W"): -2, ("Q", "Y"): -5, ("Q", "X"): 0,
        ("R", "A"): -4, ("R", "C"): -1,  ("R", "D"): -1, ("R", "E"): -8, ("R", "F"): -5, ("R", "G"):  3, ("R", "H"):  1, ("R", "I"): -3, ("R", "K"): -9, ("R", "L"): -4, ("R", "M"): -5, ("R", "N"):  2, ("R", "P"):  0, ("R", "Q"): -2, ("R", "R"):  6, ("R", "S"):  2, ("R", "T"):  0, ("R", "V"): -1, ("R", "W"): -6, ("R", "Y"): -3, ("R", "X"): 0,
        ("S", "A"): -7, ("S", "C"):  0,  ("S", "D"): -4, ("S", "E"):-10, ("S", "F"): -8, ("S", "G"): -3, ("S", "H"): -3, ("S", "I"): -9, ("S", "K"):-14, ("S", "L"): -4, ("S", "M"):-10, ("S", "N"):  3, ("S", "P"): -2, ("S", "Q"): -4, ("S", "R"):  2, ("S", "S"):  6, ("S", "T"): -6, ("S", "V"):  0, ("S", "W"):-11, ("S", "Y"): -9, ("S", "X"): 0,
        ("T", "A"): -5, ("T", "C"): -8,  ("T", "D"): -2, ("T", "E"):-10, ("T", "F"): -3, ("T", "G"):  4, ("T", "H"): -1, ("T", "I"):  6, ("T", "K"): -5, ("T", "L"): -9, ("T", "M"):  3, ("T", "N"): -4, ("T", "P"): -4, ("T", "Q"): -5, ("T", "R"):  0, ("T", "S"): -6, ("T", "T"):  8, ("T", "V"): -9, ("T", "W"): -5, ("T", "Y"): -5, ("T", "X"): 0,
        ("V", "A"): -6, ("V", "C"):  1,  ("V", "D"): -3, ("V", "E"):-13, ("V", "F"): -7, ("V", "G"): -6, ("V", "H"): -5, ("V", "I"):-12, ("V", "K"):-15, ("V", "L"):  0, ("V", "M"):-16, ("V", "N"):  0, ("V", "P"):  0, ("V", "Q"): -1, ("V", "R"): -1, ("V", "S"):  0, ("V", "T"): -9, ("V", "V"):  3, ("V", "W"):-10, ("V", "Y"):-11, ("V", "X"): 0,
        ("W", "A"):  0, ("W", "C"): -7,  ("W", "D"): -2, ("W", "E"): -6, ("W", "F"):  4, ("W", "G"): -4, ("W", "H"): -5, ("W", "I"): -5, ("W", "K"):  5, ("W", "L"): -8, ("W", "M"): -6, ("W", "N"): -8, ("W", "P"): -4, ("W", "Q"): -2, ("W", "R"): -6, ("W", "S"):-11, ("W", "T"): -5, ("W", "V"):-10, ("W", "W"):  8, ("W", "Y"): -6, ("W", "X"): 0,
        ("Y", "A"): -2, ("Y", "C"): -9,  ("Y", "D"): -2, ("Y", "E"): -3, ("Y", "F"): -4, ("Y", "G"): -2, ("Y", "H"):  3, ("Y", "I"): -8, ("Y", "K"): -8, ("Y", "L"): -9, ("Y", "M"): -9, ("Y", "N"): -5, ("Y", "P"): -5, ("Y", "Q"): -5, ("Y", "R"): -3, ("Y", "S"): -9, ("Y", "T"): -5, ("Y", "V"):-11, ("Y", "W"): -6, ("Y", "Y"):  9, ("Y", "X"): 0,
        ("X", "A"):  0, ("X", "C"):  0,  ("X", "D"):  0, ("X", "E"):  0, ("X", "F"):  0, ("X", "G"):  0, ("X", "H"):  0, ("X", "I"):  0, ("X", "K"):  0, ("X", "L"):  0, ("X", "M"):  0, ("X", "N"):  0, ("X", "P"):  0, ("X", "Q"):  0, ("X", "R"):  0, ("X", "S"):  0, ("X", "T"):  0, ("X", "V"):  0, ("X", "W"):  0, ("X", "Y"):  0, ("X", "X"): 0
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

def calculate_blosum_score(seq1, seq2, matrix, gap_open_penalty=GAP_OPEN_PENALTY, gap_penalty=GAP_PENALTY, 
                           GAP_CHAR='-', UNKNOWN_CHAR='X'):
    """
    Calculates the raw BLOSUM-like score between two aligned sequences
    using an affine gap penalty model (gap open and gap extension penalties).

    Args:
        seq1 (str): The first aligned sequence.
        seq2 (str): The second aligned sequence.
        matrix (dict): The substitution matrix (e.g., BLOSUM62 from previous context).
                       Assumed to be a dict with tuple keys, e.g., matrix[('A', 'R')] = score.
                       The function will also check for symmetric keys, e.g., matrix[('R', 'A')].
                       If the original code relied on an external `get_blosum_score(r1, r2, matrix)` 
                       helper, that call should be used instead of direct matrix access here.
        gap_open_penalty (int or float): The penalty for opening a gap (e.g., -10).
        gap_penalty (int or float): The penalty for extending a gap by one position 
                                     (referred to as gap extension penalty, e.g., -1).
        GAP_CHAR (str, optional): Character representing a gap. Defaults to '-'.
        UNKNOWN_CHAR (str, optional): Character representing an unknown residue. Defaults to 'X'.

    Returns:
        int or float: The calculated alignment score.
        
    Raises:
        ValueError: If sequences do not have the same length.
    """
    if len(seq1) != len(seq2):
        raise ValueError("Sequences must have the same length for aligned score calculation.")

    score = 0
    alignment_length = len(seq1)
    
    # State variables: True if the respective sequence is currently in a gap tract
    in_gap_s1 = False  # True if seq1 is currently in a gap that started to its left or at current position
    in_gap_s2 = False  # True if seq2 is currently in a gap that started to its left or at current position

    for i in range(alignment_length):
        res1 = seq1[i]
        res2 = seq2[i]

        # Determine if this position involves a gap character for seq1 or seq2
        is_res1_gap = (res1 == GAP_CHAR)
        is_res2_gap = (res2 == GAP_CHAR)

        if is_res1_gap and is_res2_gap:
            # Case 1: Both are gaps (e.g., '-' vs '-')
            # The original function's behavior for gap-gap was a score of 0 for the pair.
            # This position contributes 0 to the total score.
            # Both sequences are considered to be in a gapped state.
            in_gap_s1 = True
            in_gap_s2 = True
        
        elif is_res1_gap:
            # Case 2: Gap in seq1, residue in seq2 (e.g., '-' vs 'A')
            if not in_gap_s1:  # This is the start of a new gap in seq1
                score += gap_open_penalty
            score += gap_penalty  # Apply gap extension penalty for this position in seq1's gap
            
            in_gap_s1 = True   # seq1 is in a gap
            in_gap_s2 = False  # seq2 has a residue, so it's not in a gap state here
        
        elif is_res2_gap:
            # Case 3: Residue in seq1, Gap in seq2 (e.g., 'A' vs '-')
            if not in_gap_s2:  # This is the start of a new gap in seq2
                score += gap_open_penalty
            score += gap_penalty  # Apply gap extension penalty for this position in seq2's gap
            
            in_gap_s2 = True   # seq2 is in a gap
            in_gap_s1 = False  # seq1 has a residue, so it's not in a gap state here

        elif res1 == UNKNOWN_CHAR or res2 == UNKNOWN_CHAR:
            # Case 4: At least one residue is 'X' (e.g., 'X' vs 'A', 'A' vs 'X', 'X' vs 'X')
            # (and it's not a gap vs residue case, which is handled above)
            # The original function scored 0 for pairs involving 'X'.
            # This position contributes 0 to the total score.
            in_gap_s1 = False # Not in a gap state
            in_gap_s2 = False # Not in a gap state
        
        else:
            # Case 5: Standard residue pair (e.g., 'A' vs 'C')
            pair_value = 0  # Default if pair is not found in the matrix
            
            # Attempt to get score from the matrix (assuming tuple-keyed dictionary)
            if (res1, res2) in matrix:
                pair_value = matrix[(res1, res2)]
            elif (res2, res1) in matrix: # Check for symmetric key (e.g., if matrix stores ('C','A') but not ('A','C'))
                pair_value = matrix[(res2, res1)]
            # else:
                # Optional: handle or log missing pairs if matrix might be incomplete.
                # For example: print(f"Warning: Score for pair ({res1}, {res2}) not found. Using {pair_value}.")
            
            score += pair_value
            
            # Since both are residues, neither sequence is in a gap state at this point
            in_gap_s1 = False 
            in_gap_s2 = False
            
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
                        f_out.write(f"{len(cluster)} " + " ".join(cluster) + "\n")
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
