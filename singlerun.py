#!/usr/bin/python
import argparse, math, os, sys
from alignment import Alignment
from params import parse_params
from scorers.base import get_scorer
from utils.bio import get_column


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
            scores = [None] * len(alignment.msa[0])
        all_scores.append(scores)
    # Rearrange to list of scores, grouped by site
    score_tups = zip(*all_scores)
    return score_tups


################################################################################
# Input/output
################################################################################

def prepare_header(scorers):
    header = []
    for scorer in scorers:
        header.append("# %s" % scorer.name)
        params = scorer.get_params()
        for k, v in params:
            header.append("# \t%s: %s" % (k, v))
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
        f.write("%d\t%s\t%s\n" % (i+1, site, "\t".join(map(_write_score_helper, score_tup))))


def read_scores(fname):
    prevline = None
    score_tups = []
    with open(fname) as f:
        for line in f:
            if line.startswith('#'):
                prevline = line
                continue
            line = line.strip()
            if not line:
                continue
            if prevline:
                scorer_names = prevline.strip().split()[2:]
            prevline = None
            score_tups.append(map(_read_score_helper, line.split()[2:]))
    return score_tups


def _write_score_helper(score):
    if isinstance(score, float):
        return str(round(score,4))
    elif score is None:
        return '-'
    return str(score)

def _read_score_helper(field):
    if field == '-':
        return None
    return float(field)


################################################################################
# Cmd line driver
################################################################################

def parse_args():
    parser = argparse.ArgumentParser(
        description="Produce conservation scores for an alignment file.")

    parser.add_argument('scorer_name',
        help="conservation estimation method")
    parser.add_argument('align_file', 
        help="path to alignment file to score")

    parser.add_argument('-a', dest='align_params', action='append', default=[],
        help="parameters associcated with align_file, can specify multiple. Specify as '-a inputName=inputValue', e.g. '-a tree_file=tree.txt'")
    parser.add_argument('-p', dest='scorer_params', action='append', default=[],
        help="parameters to pass to the scorer, can specify multiple. Specify as '-p paramName=paramValue', e.g. '-p use_gap_penalty=1")

    parser.add_argument('-o', dest='out_fname', type=str, 
        help="name of output file. Default prints to screen")
    parser.add_argument('-r', dest='reference_seq', type=str,
        help="reference sequence to print scores in reference to (ignoring gaps). Default prints the entire column.")

    args = parser.parse_args()
    args.align_params = parse_params(args.align_params)
    args.scorer_params = parse_params(args.scorer_params)

    return args


def main():
    args = parse_args()

    # Get scorer(s)
    scorers = [get_scorer(args.scorer_name, **args.scorer_params)]

    # Get alignment and supplementary inputs
    alignment = Alignment(args.align_file, **args.align_params)

    score_tups = compute_scores(alignment, scorers)

    header = prepare_header(scorers)

    scorer_names = [args.scorer_name]
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