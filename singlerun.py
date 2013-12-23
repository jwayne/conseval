#!/usr/bin/env python

import argparse, math, os, sys
from alignment import Alignment
from scorer import get_scorer
from utils.bio import get_column


################################################################################
# Input/output
################################################################################

DEFAULT_SCORER = 'mayrose04'
def parse_scorer_names(scorer_names):
    if not scorer_names:
        return [DEFAULT_SCORER]
    return scorer_names


def compute_scores(alignment, scorers, all_scores=None):
    # Compute list of scores, grouped by scorer
    if not all_scores:
        all_scores = []
    for scorer in scorers:
        try:
            scores = scorer.score(alignment)
        except Exception, e:
            import traceback
            sys.stderr.write("Error scoring %s via %s\n" %
                (alignment.align_file, type(scorer).__name__))
            traceback.print_exc()
            scores = [None] * len(alignment.msa)
        all_scores.append(scores)
    # Rearrange to list of scores, grouped by site
    score_tups = zip(*all_scores)
    return score_tups


def parse_score(score):
    if isinstance(score, float):
        return str(round(score,4))
    elif score is None:
        return '-'
    return str(score)


def write_scores(alignment, score_tups, scorer_names, f=sys.stdout,
        header=None, reference_seq=None, includes_testset=False):

    if header:
        # Check if header lines start with '#'?
        f.write(header)
        f.write("\n\n")

    # handle printing of aa's from 1 sequence, rather than aa's from entire alignment
    ref_seq_num = None
    if reference_seq:
        if reference_seq in alignment.names:
            ref_seq_num = alignment.names.index(reference_seq)
        else:
            sys.stderr.write("Reference sequence '%s' not found in alignment. Using default output format\n" % reference_seq)
    if ref_seq_num is not None:
        f.write("# reference sequence: %s\n" % reference_seq)

    if includes_testset:
        if not len(score_tups[0]) == len(scorer_names)+1:
            raise ValueError("Mismatch between inputs 'score_tups' and 'scorer_names'")
        f.write("# col\tsite\ttestset\t%s\n" % "\t".join(scorer_names))
    else:
        if not len(score_tups[0]) == len(scorer_names):
            raise ValueError("Mismatch between inputs 'score_tups' and 'scorer_names'")
        f.write("# col\tsite\t%s\n" % "\t".join(scorer_names))

    # print scores
    for i, score_tup in enumerate(score_tups):
        if ref_seq_num is not None:
            cur_aa = get_column(i, alignment.msa)[ref_seq_num]
            if cur_aa == '-': continue
            site = cur_aa
        else:
            site = "".join(get_column(i, alignment.msa))
        f.write("%d\t%s\t%s\n" % (i+1, site, "\t".join(map(parse_score,score_tup))))



################################################################################
# Cmd line driver
################################################################################

def parse_args():
    parser = argparse.ArgumentParser(description="Produce conservation scores for an alignment file.")

    parser.add_argument('align_file')
    parser.add_argument('-a', dest='align_params', action='append', default=[],
        help="parameters associcated with align_file, can specify multiple. Specify as '-a inputName=inputValue', e.g. '-a tree_file=tree.txt'")

    # TODO: accept a single scorer only? leave running multiple scorers to experiment.py?
    parser.add_argument('-s', dest='scorer_names', action='append', default=[],
        help="conservation estimation method, can specify multple. Default='%s'" % DEFAULT_SCORER)
    parser.add_argument('-p', dest='params', action='append', default=[],
        help="parameters to pass to the scorer, can specify multiple. Specify as '-p paramName=paramValue', e.g. '-p gap_penalty=1")

    parser.add_argument('-o', dest='out_fname', type=str, 
        help="name of output file. Default prints to screen")
    parser.add_argument('-r', dest='reference_seq', type=str,
        help="reference sequence to print scores in reference to (ignoring gaps). Default prints the entire column.")

    args = parser.parse_args()
    args.align_params = dict(val.split("=") for val in args.align_params)
    args.params = dict(val.split("=") for val in args.params)
    return args


def main():
    args = parse_args()
    align_file = args.align_file

    # Get scorer(s)
    scorer_names = parse_scorer_names(args.scorer_names)
    scorers = [get_scorer(name, **args.params) for name in scorer_names]

    # Get alignment and supplementary inputs
    alignment = Alignment(align_file, **args.align_params)

    score_tups = compute_scores(alignment, scorers)

    header = []
    for scorer in scorers:
        header.append("# %s params:" % scorer.name)
        params = scorer.get_params()
        for k in sorted(params.keys()):
            header.append("# \t%s: %s" % (k, params[k]))
    header = "\n".join(header)

    if args.out_fname:
        with open(args.out_fname, 'w'):
            write_scores(alignment, score_tups, scorer_names, f,
                header=header,
                reference_seq=args.reference_seq)
    else:
        write_scores(alignment, score_tups, scorer_names,
            header=header,
            reference_seq=args.reference_seq)


if __name__ == "__main__":
    main()
