import os
from Bio import Phylo, SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


#####
# Conversion tools
#####

def get_phylotree(alignment, n_bootstrap=0, overwrite=False):
    """
    Get the phylo tree corresponding to `alignment`.  If no tree, compute one
    and cache to disk.
    """
    fname_phy = '.'.join(alignment.align_file.split('.')[:-1]) + '.phy'
    fname_tree = fname_phy + '_phyml_tree.txt'
    tree = None
    if not overwrite and os.path.exists(fname_tree) and os.path.getsize(fname_tree):
        tree = read_phylotree(fname_tree)
        # Check that the cached tree matches the alignment. If not, re-compute.
        tree_names = set(clade.name for clade in tree.get_terminals())
        msa_names = set(alignment.names)
        if tree_names != msa_names:
            tree = None
    if not tree:
        tree = _compute_phylotree(alignment, fname_phy, fname_tree, n_bootstrap)
    return tree


def read_phylotree(fname_tree):
    return Phylo.read(fname_tree, "newick")


def _compute_phylotree(alignment, fname_phy, fname_tree, n_bootstrap):
    """
    Use PhyML to compute the tree for fname_aln.
    Does not compute tree if treefile exists and overwrite=False (default).
    Return the filename of the tree.
    """
    records = []
    for row,name in zip(alignment.msa, alignment.names):
        records.append(SeqRecord(Seq(row), id=name, description=name))
    with open(fname_phy, "w") as f_out:
        SeqIO.write(records, f_out, "phylip-relaxed")
    os.system("phyml -i %s -d aa -b %d --quiet --no_memory_check" % (fname_phy, n_bootstrap))
    return read_phylotree(fname_tree)


#####

if __name__ == "__main__":
    import sys
    from alignment import Alignment
    print get_phylotree(Alignment(sys.argv[1]))
