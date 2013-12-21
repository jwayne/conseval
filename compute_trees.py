#!/usr/bin/python
import glob
from multiprocessing import Pool, TimeoutError
import os
import sys
import time

from alignment import Alignment
from utils_parallel import parallelize


def compute_tree(fn):
    t = time.time()
    aln = Alignment(fn)
    aln.get_phylotree()
    return time.time() - t


def compute_trees(fns, out_file):
    it = parallelize(compute_tree, fns, timeout=1800)

    f = open(out_file, 'w')
    for arg, result in it:
        if result is None:
            f.write("Timeout\t%s\n" % arg)
            f.flush()
        else:
            f.write("%7.2f\t%s\n" % (result, arg))
            f.flush()
    f.close()


if __name__ == "__main__":
    fns = glob.glob(os.path.join(sys.argv[1], "*.aln"))
    out_file = sys.argv[2]
    compute_trees(fns, out_file)
