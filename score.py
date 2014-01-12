#!/usr/bin/python
"""
Score a single alignment using a single scorer.  Print scores in
a pleasant human-readable format.
"""
import argparse
import sys
from conseval.alignment import Alignment
from conseval.io import write_score_helper, read_score_helper, list_scorer_params, parse_params
from conseval.scorer import get_scorer, get_scorer_cls
from conseval.utils.bio import get_column
from conseval.utils.general import get_all_module_names



################################################################################
# Input/output
################################################################################

def write_scores(alignment, score_tups, scorer_names, f=sys.stdout, header=None):
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

    # print scores
    f.write("# i\tcolumn\t%s\n" % "\t".join(scorer_names))
    for i, score_tup in enumerate(score_tups):
        site = "".join(get_column(i, alignment.msa))
        f.write("%d\t%s\t%s\n" % (i+1, site, "\t".join(map(write_score_helper, score_tup))))


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
            score_tups.append(map(read_score_helper, line.split()[2:]))
    return score_tups


def list_alignment_paramdefs():
    print "================"
    print "Alignment params:"
    print "================"
    print "\n".join("    %s" % pd for pd in Alignment.params.param_defs)
    print ""


def list_scorer_paramdefs(scorer_names=None):
    if not scorer_names:
        scorer_names = get_all_module_names('scorers')
    for scorer_name in scorer_names:
        scorer_cls = get_scorer_cls(scorer_name)
        x = "Scorer params for %s:" % scorer_name
        print ("="*len(x))
        print x
        print ("="*len(x))
        print "\n".join("    %s" % pd for pd in scorer_cls.params.param_defs)
        print ""


################################################################################
# Cmd line driver
################################################################################

def parse_args():
    parser = argparse.ArgumentParser(
        description="Score the conservation of a single alignment file, using a single scorer.",
        usage="%(prog)s [-h] [-l] scorer_name align_file [-a ALIGN_PARAMS] [-p SCORER_PARAMS]")

    parser.add_argument('-l', dest='list_params', nargs=argparse.REMAINDER,
        help="list parameters for all scorers (score.py -l), or for certain defined scorers (score.py -l intrepid rate4site_eb)")

    parser.add_argument('scorer_name', nargs="?",
        help="conservation estimation method")
    parser.add_argument('align_file', nargs="?",
        help="path to alignment file to score")

    parser.add_argument('-a', dest='align_params', action='append', default=[],
        help="parameters associcated with align_file, can specify multiple. Specify as '-a inputName=inputValue', e.g. '-a tree_file=tree.txt'")
    parser.add_argument('-p', dest='scorer_params', action='append', default=[],
        help="parameters to pass to the scorer, can specify multiple. Specify as '-p paramName=paramValue', e.g. '-p use_gap_penalty=1")

    args = parser.parse_args()
    if args.list_params:
        list_alignment_paramdefs()
        list_scorer_paramdefs(args.list_params)
        sys.exit(0)
    elif not args.scorer_name or not args.align_file:
        parser.print_usage()
        sys.stderr.write("%s: error: too few arguments\n" % sys.argv[0])
        sys.exit(1)
    args.align_params = parse_params(args.align_params)
    args.scorer_params = parse_params(args.scorer_params)

    return args


def main():
    args = parse_args()

    # Get scorer
    scorer = get_scorer(args.scorer_name, **args.scorer_params)

    # Get alignment and supplementary inputs
    alignment = Alignment(args.align_file, **args.align_params)

    # Score
    scores = scorer.score(alignment)

    # Output
    header = list_scorer_params(scorer)
    write_scores(alignment, zip(scores), [args.scorer_name], header=header, f=sys.stdout)


if __name__ == "__main__":
    main()
