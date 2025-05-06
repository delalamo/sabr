import os
import subprocess
import sys

from Bio import PDB

THREE_TO_ONE_MAP = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V"
}

parser = PDB.PDBParser(QUIET=True)

structdir = "/home/delalamo/sabdab_structures/parsed/pdb_segments/"
foldseekdir = "/home/delalamo/sabdab_structures/parsed/foldseek/"
subdirs = ["camelid_heavy", "human_heavy", "human_light", "mouse_heavy", "mouse_light"]

problem_seqs = ["1jto", "1jtp", "1jtt", "1mel", "1sjv"]

for subdir in subdirs:
    # structs = {}
    # for i, file in enumerate(os.listdir(os.path.join(structdir, subdir))):
    #     name = os.path.join(structdir, subdir, file)
    #     print(i, name)
    #     structs[name] = parser.get_structure("TEMP", name)
        
    struct_file = os.path.join(foldseekdir, subdir, f"{subdir}_ss")
    name_file = os.path.join(foldseekdir, subdir, f"{subdir}_h")
    seq_file = os.path.join(foldseekdir, subdir, f"{subdir}")


    all_structs, all_names, all_seqs = [], [], []
    with open(struct_file, 'rb') as f:
        for line in f:
            line = line.rstrip(b'\x00').lstrip(b'\x00')
            if len(line.rstrip(b'\x00')) > 1:
                all_structs.append(line.decode('utf-8').strip())
    with open(name_file, 'rb') as f:
        for line in f:
            line = line.rstrip(b'\x00').lstrip(b'\x00')
            if len(line.strip()) > 1:
                all_names.append(line.decode('utf-8').strip())
    with open(seq_file, 'rb') as f:
        for line in f:
            line = line.rstrip(b'\x00').lstrip(b'\x00')
            if len(line.strip()) > 1:
                all_seqs.append(line.decode('utf-8').strip())

    assert len(all_seqs) == len(all_names)
    assert len(all_seqs) == len(all_structs)

    print(len(all_structs))
    print(len(all_seqs))
    print(len(all_names))

    mismatches_seqs = []

    seqs = []

    start = 5
    cutoff = 128
    if "light" in subdir:
        cutoff = 127
    for j, struct_tokens, name, seq in zip(range(len(all_names))[start:], all_structs[start:], all_names[start:], all_seqs[start:]):
        if name[:4] in problem_seqs:
            continue
        filename = os.path.join(structdir, subdir, f"{name}.pdb")
        struct = parser.get_structure("TEMP", filename)
        # print(filename)
        # if filename not in structs:
        #     print("File not found: ", filename)
        #     missed_files.append(filename)
        #     continue
        # struct = structs[filename]
        aln = ["-"] * 128
        realseq = "".join([THREE_TO_ONE_MAP[r.get_resname()] for r in struct.get_residues()])
        if realseq != seq:
            mismatches_seqs.append(filename)
            continue

        idxs_real = [r.get_id()[1:] for r in struct.get_residues()]
        idxs_test = []
        cmd = ["ANARCI", "--sequence", struct_tokens]
        out = subprocess.run(cmd, capture_output=True, text=True, check=True)
        species = None
        chain = None
        startpos = None
        for i, line in enumerate(out.stdout.split("\n")):
            if i == 5:
                species = line.split("|")[1]
                chain = line.split("|")[2]
                startpos = int(line.split("|")[5])
            if not line.startswith("#"):
                splitline = line.split()
                if len(splitline) < 3 or line[-1] == "-":
                    continue
                res_id = (int(splitline[1]), ' ')
                if len(splitline) == 4:
                    res_id = (int(splitline[1]), splitline[2])
                idxs_test.append(res_id)
        mismatch = False
        for real, test in zip(idxs_real[startpos:], idxs_test):
            if real != test:
                mismatch = True
                break
        mismatch_species = False
        real_species = subdir.split("_")[0]
        real_chain = subdir[-5].upper()
        if species != real_species:
            mismatch_species = True
        mismatch_chain = False
        if chain != real_chain:
            mismatch_chain = True
        if mismatch:
            print(f"mismatch found: {name} {species} {chain} ({real_species} {real_chain}) {len(idxs_real)} vs {len(idxs_test)} {struct_tokens}")
            for real, test in zip(idxs_real, idxs_test):
                print(real, test)
            sys.exit()
print("Skipping mismatched seqs")
for mismatch in mismatches_seqs:
    print(mismatch)

# problem PDBs
# 1jto