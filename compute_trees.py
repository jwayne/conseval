#!/usr/bin/python
import argparse
import glob
from multiprocessing import Pool, TimeoutError
import os
import sys
import time

from alignment import Alignment
from dataset_config import DATASET_CONFIGS
from utils import parallelize
from utils.general import get_timestamp


def compute_tree(fn):
    t = time.time()
    aln = Alignment(fn)
    aln.get_phylotree()
    return time.time() - t


def compute_trees(fns, out_file):
    it = parallelize.imap_unordered(compute_tree, fns, timeout=1800)

    with open(out_file, 'a') as f:
        try:
            f.write("# Start time: %s\n" % get_timestamp())
            f.write("# Total num alignments: %d\n" % len(fns))
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
    parser.add_argument('-f', dest="force", action="store_true",
        help="overwrite tree files if not matching. Default doesn't overwrite if tree file exists.")
    args = parser.parse_args()
    
    dataset_config = DATASET_CONFIGS[args.dataset_name]
    fns = dataset_config.get_filenames()
    aln_dir = dataset_config.aln_dir

    if not args.force:
        fns_old = fns
        fns = []
        for fn in fns_old:
            align_file = os.path.join(aln_dir, fn)
            # TODO: not just aln
            tree_file = align_file[:-3] + "phy_phyml_tree.txt"
            if not os.path.exists(tree_file):
                fns.append(align_file)
    print "# Total num alignments: %d\n" % len(fns)
    compute_trees(fns, args.out_file)
