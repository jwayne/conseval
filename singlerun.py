#!/usr/bin/env python

import argparse, math, os, sys
from alignment import Alignment
from params import parse_params
from scorer import get_scorer, get_all_scorer_names, get_scorer_cls
from utils.bio import get_column


################################################################################
# Input/output
################################################################################

def compute_scores(alignment, scorers):
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


def prepare_header(scorers):
    header = []
    for scorer in scorers:
        header.append("# %s" % scorer.name)
        params = scorer.get_params()
        for k in sorted(params.keys()):
            header.append("# \t%s: %s" % (k, params[k]))
    return "\n".join(header) + "\n"


def write_scores(alignment, score_tups, scorer_names, f=sys.stdout,
        header=None, reference_seq=None):

    # sanity check
    if not len(score_tups[0]) == len(scorer_names):
        raise ValueError("Mismatch between inputs 'score_tups' and 'scorer_names'")

    f.write("# Alignment: %s\n" % alignment.align_file)
    f.write("# Num sites: %d\n" % len(alignment.msa[0]))
    f.write("# Num sequences: %d\n" % len(alignment.names))
    if alignment.filtered:
        f.write("# Num sequences before filtering: %d\n" % alignment.orig_num_sequences)
    f.write("\n")

    if header:
        # Check if header lines start with '#'?
        f.write(header)
        f.write("\n")

    # handle printing of aa's from 1 sequence, rather than aa's from entire alignment
    ref_seq_num = None
    if reference_seq:
        if reference_seq in alignment.names:
            ref_seq_num = alignment.names.index(reference_seq)
        else:
            sys.stderr.write("Reference sequence '%s' not found in alignment. Using default output format\n" % reference_seq)
    if ref_seq_num is not None:
        f.write("# reference sequence: %s\n" % reference_seq)

    # print scores
    f.write("# col\tsite\t%s\n" % "\t".join(scorer_names))
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
    parser = argparse.ArgumentParser(
        description="Produce conservation scores for an alignment file.",
        usage="%(prog)s scorer_name align_file [-h] [-l] [options]")

    parser.add_argument('scorer_name', nargs='?',
        help="conservation estimation method")
    parser.add_argument('align_file', nargs='?',
        help="path to alignment file to score")

    parser.add_argument('-l', '--list', dest='list_params', action='store_true',
        help="show list of available params for the -a, -p flags")

    parser.add_argument('-a', dest='align_params', action='append', default=[],
        help="parameters associcated with align_file, can specify multiple. Specify as '-a inputName=inputValue', e.g. '-a tree_file=tree.txt'")
    parser.add_argument('-p', dest='scorer_params', action='append', default=[],
        help="parameters to pass to the scorer, can specify multiple. Specify as '-p paramName=paramValue', e.g. '-p gap_penalty=1")

    parser.add_argument('-o', dest='out_fname', type=str, 
        help="name of output file. Default prints to screen")
    parser.add_argument('-r', dest='reference_seq', type=str,
        help="reference sequence to print scores in reference to (ignoring gaps). Default prints the entire column.")

    args = parser.parse_args()
    args.align_params = parse_params(args.align_params)
    args.scorer_params = parse_params(args.scorer_params)

    if args.list_params:
        print "Alignment params:"
        print "\n".join("\t%s" % pd for pd in Alignment.PARAMS.param_defs)
        print ""
        if args.scorer_name:
            scorer_names = [args.scorer_name]
        else:
            scorer_names = get_all_scorer_names()
        for scorer_name in scorer_names:
            print "Scorer params for %s:" % scorer_name
            scorer_cls = get_scorer_cls(scorer_name)
            print "\n".join("\t%s" % pd for pd in scorer_cls.PARAMS.param_defs)
            print ""
        sys.exit(0)

    if not args.scorer_name or not args.align_file:
        parser.print_help()
        sys.exit(0)

    return args


def main():
    args = parse_args()

    # Get scorer(s)
    scorers = [get_scorer(args.scorer_name, **args.scorer_params)]

    # Get alignment and supplementary inputs
    alignment = Alignment(args.align_file, **args.align_params)

    score_tups = compute_scores(alignment, scorers)

    header = prepare_header(scorers)

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
