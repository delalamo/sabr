#!/usr/bin/env python3
"""Constants and configuration values for SAbR.

This module defines constants used throughout the SAbR package including:
- Neural network embedding dimensions
- IMGT numbering scheme definitions
- Amino acid mappings
- Alignment parameters
"""

##################
# SoftAlign and MPNN constants

# SoftAlign params
EMBED_DIM = 64
N_MPNN_LAYERS = 3

# Backbone CB reconstruction params
CB_BOND_LENGTH = 1.522
CB_BOND_ANGLE = 1.927
CB_DIHEDRAL = -2.143
PEPTIDE_BOND_LENGTH = 1.33

# Cutoff for gap detection
PEPTIDE_BOND_MAX_DISTANCE = 2 * PEPTIDE_BOND_LENGTH
MAX_SELECTED_RESIDUES = 1024

# Backbone indices for MPNN
BACKBONE_N_IDX = 0
BACKBONE_CA_IDX = 1
BACKBONE_C_IDX = 2
BACKBONE_CB_IDX = 3

# Gap scores (from NPZ)
SW_GAP_EXTEND = -0.175027
SW_GAP_OPEN = -2.525591

# SoftAlign params
DEFAULT_TEMPERATURE = 1e-4
CHAIN_TYPES = ("H", "K", "L")
SCFV_CHAIN_TYPES = ("H:K", "H:L", "K:H", "L:H")
NOISE_LEVELS = (0.0, 0.2, 0.5, 1.0, 2.0)
NUMBERING_SCHEMES = ("imgt", "chothia", "kabat", "martin", "aho", "wolfguy")

##################
# IMGT constants

IMGT_MAX_POSITION = 128

# DE loop positions 81-84
FR3_POS81_COL = 80
FR3_POS82_COL = 81
FR3_POS83_COL = 82
FR3_POS84_COL = 83

# C-terminus correction position
C_TERMINUS_ANCHOR_POSITION = 124

# CDR loop definitions (inclusive)
IMGT_LOOPS = {
    "CDR1": (27, 38),
    "CDR2": (56, 65),
    "CDR3": (105, 117),
}

# Anchor points for deterministic renumbering
# These are spatially conserved residue with high accuracy
# Not necessarily CDR termini!
CDR_ANCHORS = {
    "CDR1": (23, 40),
    "CDR2": (54, 67),
    "CDR3": (104, 118),
}

##################
# Misc

AA_3TO1 = {
    "ALA": "A",
    "CYS": "C",
    "ASP": "D",
    "GLU": "E",
    "PHE": "F",
    "GLY": "G",
    "HIS": "H",
    "ILE": "I",
    "LYS": "K",
    "LEU": "L",
    "MET": "M",
    "ASN": "N",
    "PRO": "P",
    "GLN": "Q",
    "ARG": "R",
    "SER": "S",
    "THR": "T",
    "VAL": "V",
    "TRP": "W",
    "TYR": "Y",
}
