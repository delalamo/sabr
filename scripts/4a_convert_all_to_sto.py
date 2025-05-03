import os

from Bio import SeqIO


def write_stockholm(seqs, ID, outfile):
    print("# STOCKHOLM 1.0", file=outfile)
    print(f"#=GF ID {ID}", file=outfile)

    pad_length = max(list(map(len, seqs)))
    for i, s in enumerate(seqs):
        print(s.replace(" ", "_").ljust(pad_length), seqs[i].replace(".", "-"), file=outfile)
    print("#=GC RF".ljust(pad_length), "x"*len(seqs[0]), file=outfile)
    print("//", file=outfile)

def output_stockholm_all(all_seqs, path):
    filename = os.path.join(path, "ALL.sto")
    with open(filename, "w") as f:
        for (species, chain_type), seqs in all_seqs.items():
            write_stockholm(seqs, f"{species}_{chain_type}", f)

def get_seqs_from_a3m(file):
    seqs = []
    for record in SeqIO.parse(file, "fasta"):
        seq = str(record.seq)
        if seq not in seqs:
            seqs.append(str(record.seq))
    return seqs

def main():
    allseqs = {}
    for file in os.listdir("/home/delalamo/sabr/a3ms/"):
        if not file.endswith("a3m"):
            continue
        species, chain = file.split("_")
        seqs = get_seqs_from_a3m(os.path.join("/home/delalamo/sabr/a3ms/", file))
        allseqs[(species, chain[0])] = seqs
    output_stockholm_all(allseqs, "/home/delalamo/sabr/stos/")

if __name__ == "__main__":
    main()
