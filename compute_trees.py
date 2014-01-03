#!/usr/bin/python
import argparse
import glob
from multiprocessing import Pool, TimeoutError
import os
import sys
import time

from alignment import Alignment
from dataset_config import DATASET_CONFIGS
from phylotree import read_phylotree, check_phylotree
from utils import parallelize
from utils.general import get_timestamp


def compute_tree(align_file):
    t = time.time()
    aln = Alignment(align_file)
    aln.get_phylotree()
    return time.time() - t


def compute_trees(align_files, out_file):
    it = parallelize.imap_unordered(compute_tree, align_files, timeout=1800)

    with open(out_file, 'a') as f:
        try:
            f.write("# Start time: %s\n" % get_timestamp())
            f.write("# Total num alignments: %d\n" % len(align_files))
            count = 0
            for arg, result in it:
                if result is None:
                    timestr = "Timeout"
                else:
                    timestr = "%7.2f" % result
                f.write("%d\t%s\t%s\n" % (count, result, arg))
                f.flush()
                count += 1
        finally:
            f.write("# End time: %s\n" % get_timestamp())


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Compute phylogenetic trees on a set of files.")
    parser.add_argument("dataset_name")
    parser.add_argument("out_file")
    parser.add_argument('--check-tree-contents', action="store_true",
        help="check if contents of tree files match the alignment, and re-compute trees if they don't. Default only checks if the corresponding tree file for the alignment exists.")
    args = parser.parse_args()
    
    dataset_config = DATASET_CONFIGS[args.dataset_name]
    align_files = dataset_config.get_align_files()

    if not args.check_tree_contents:
        align_files_old = align_files
        align_files = []
        for align_file in align_files_old:
            # TODO: handle extensions besides .aln
            tree_file = align_file[:-3] + "phy_phyml_tree.txt"
            if not os.path.exists(tree_file):
                align_files.append(align_file)
            else:
                alignment = Alignment(align_file)
                try:
                    tree = read_phylotree(tree_file)
                except:
                    align_files.append(align_file)
                    continue
                if not check_phylotree(alignment, tree):
                    align_files.append(align_file)
    print "# Total num alignments: %d\n" % len(align_files)
    compute_trees(align_files, args.out_file)
