import os
import sys
from Bio import Phylo, SeqIO


def convert_fname_aln2phy(fname_aln):
    return fname_aln[:-3] + "phy"

def convert_file_aln2phy(fname_aln):
    """
    Use BioPython to convert the aln file to a phy file.
    """
    fname_phy = convert_fname_aln2phy(fname_aln)
    with open(fname_aln, "rU") as f_in:
        with open(fname_phy, "w") as f_out:
            SeqIO.convert(f_in, "clustal", f_out, "phylip-relaxed")
    return fname_phy


def convert_fname_aln2tree(fname_aln):
    return fname_aln[:-3] + "phy_phyml_tree.txt"

def run_phyml(fname_phy, n_bootstrap):
    os.system("phyml -i %s -d aa -b %d" % (fname_phy, n_bootstrap))

def compute_tree(fname_aln, n_bootstrap=5, overwrite=False):
    """
    Use PhyML to compute the tree for fname_aln.
    Does not compute tree if treefile exists and overwrite=False (default).
    Return the filename of the tree.
    """
    fname_tree = convert_fname_aln2tree(fname_aln)
    if overwrite or not os.path.exists(fname_tree):
        fname_phy = convert_fname_aln2phy(fname_aln)
        if not os.path.exists(fname_phy):
            convert_file_aln2phy(fname_aln)
        run_phyml(fname_phy, n_bootstrap)
    return fname_tree

def get_tree(fname_aln, n_bootstrap=0, overwrite=False):
    fname_tree = compute_tree(fname_aln, n_bootstrap, overwrite)
    return Phylo.read(fname_tree, "newick")


if __name__ == "__main__":
    print get_tree(sys.argv[1])
