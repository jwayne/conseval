#!/usr/bin/env python

import argparse, math, os, sys
from alignment import Alignment
from scorer import get_scorer, parse_scorer_names, DEFAULT_SCORER
from utils import get_column


################################################################################
# Execution process
################################################################################

def print_scores(scores, alignment, scorer_names, kwargs,
        out_fname=None, out_seq=None):

    if out_fname:
        try:
            f = open(out_fname, 'w')
        except IOError, e:
            sys.stderr.write("Could not open %s for output. Printing results to standard out.\n" % out_fname)
    else:
        f = sys.stdout

    # Print kwargs
    for pkey in sorted(kwargs.keys()):
        f.write("# %s - %s\n" % (pkey, kwargs[pkey]))

    # handle print of output relative to specific sequence
    ref_seq_num = None
    if out_seq:
        if out_seq not in alignment.names:
            sys.stderr.write("Sequence %s not found in alignment. Using default output format...\n" % out_seq)
            out_seq = 0
        else:
            ref_seq_num = alignment.names.index(out_seq)

    if out_seq:
        f.write("# reference sequence: %s\n" % out_seq)
        f.write("# align_column_number\tamino acid\n%s\n" % "\t".join(scorer_names))
    else:
        f.write("# align_column_number\tcolumn\n%s\n" % "\t".join(scorer_names))

    # print scores
    scores = zip(scores)
    for i, score in enumerate(scores):
        if ref_seq_num is not None:
            cur_aa = get_column(i, alignment.msa)[ref_seq_num]
            if cur_aa == '-': continue
            annotation = cur_aa
        else:
            annotation = "".join(get_column(i, alignment.msa))
        f.write("%d\t%s\t%s\n" % (i+1, annotation, "\t".join(("%.5f" % s for s in score))))

    if out_fname:
        f.close()



################################################################################
# Cmd line driver
################################################################################

def parse_args():
    parser = argparse.ArgumentParser(description="Produce conservation scores for an alignment file.")

    parser.add_argument('align_file')
    parser.add_argument('-x', dest='extra_inputs', action='append', default=[],
        help="extra inputs associcated with align_file, can specify multiple. Specify as '-x inputName=inputValue', e.g. '-x tree_file=tree.txt'")

    parser.add_argument('-s', dest='scorer_names', action='append', default=[],
        help="conservation estimation method, can specify multple. Default='%s'" % DEFAULT_SCORER)
    parser.add_argument('-p', dest='params', action='append', default=[],
        help="parameters to pass to the scorer, can specify multiple. Specify as '-p paramName=paramValue', e.g. '-p gap_penalty=1")

    parser.add_argument('-o', dest='out_fname', type=str, 
        help="name of output file. Default prints to screen")
    parser.add_argument('-a', dest='out_seq', type=str,
        help="reference sequence to print scores in reference to (ignoring gaps). Default prints the entire column.")

    args = parser.parse_args()
    args.extra_inputs = dict(val.split("=") for val in args.extra_inputs)
    args.params = dict(val.split("=") for val in args.params)
    return args


def main():
    args = parse_args()
    align_file = args.align_file

    # Get scorer(s)
    scorer_names = parse_scorer_names(args.scorer_names)
    scorers = [get_scorer(name, **args.params) for name in scorer_names]

    # Get alignment and supplementary inputs
    alignment = Alignment(align_file, **args.extra_inputs)

    # Compute list of scores, grouped by scorer
    all_scores = []
    for scorer in scorers:
        try:
            scores = scorer.score(alignment)
        except Exception, e:
            import traceback
            sys.stderr.write("Error scoring %s via %s\n" %
                (alignment.align_file, type(scorer).__name__))
            traceback.print_exc()
            scores = [scorer.DEFAULT_SCORE] * len(alignment.msa)
        all_scores.append(scores)

    # Rearrange to list of scores, grouped by site
    all_scores = zip(*all_scores)

    kwargs = vars(args)
    print_scores(scores, alignment, scorer_names, kwargs,
        args.out_fname, args.out_seq)


if __name__ == "__main__":
    main()
