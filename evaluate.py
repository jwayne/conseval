#!/usr/bin/python
import argparse
import imp
import os
import sys

from conseval.alignment import Alignment
from conseval.datasets import DATASET_CONFIGS
from conseval.io import read_batchscores, parse_params, OUTPUT_DIR
from conseval.utils.bio import get_column
from conseval.utils.general import get_timestamp



def get_batchscore_dir(dataset_name):
    """
    Get the directory where batchscore.py output scores for `dataset_name`
    are stored.
    """
    return os.path.join(OUTPUT_DIR, "batchscore-%s" % dataset_name)


def get_batchscores(dataset_name, batchscore_ids=[], align_files_only=False):
    """
    Useful for evaluators.
    Get an iterator on (alignment, scores_col) where scores_col consists
    of lists of scores for each id in `batchscore_ids`.  This iterator is over
    all alignments in `dataset_name`.
    """
    # Sanity check.
    ds_dir = get_batchscore_dir(dataset_name)
    if not os.path.exists(ds_dir):
        raise IOError("%s for dataset %r does not exist"
                % (ds_dir, dataset_name))
    for batchscore_id in batchscore_ids:
        sc_dir = os.path.join(ds_dir, batchscore_id)
        if not os.path.exists(sc_dir):
            raise IOError("%s for dataset %r, scorer %r does not exist"
                    % (sc_dir, dataset_name, batchscore_id))

    dataset_config = DATASET_CONFIGS[dataset_name]
    align_files = dataset_config.get_align_files()

    # Be particular about which alignments we can evaluate.
    afs = []
    for align_file in align_files:
        alignment = Alignment(align_file)
        n_gapped_cols = 0
        for i in xrange(len(alignment.msa[0])):
            col = get_column(i, alignment.msa)
            if col.count('-') > len(col) / 2:
                n_gapped_cols += 1
        if n_gapped_cols > len(alignment.msa[0]) / 2:
            continue
        include = True
        for batchscore_id in batchscore_ids:
            sc_dir = os.path.join(ds_dir, batchscore_id)
            out_file = dataset_config.get_out_file(align_file, sc_dir)
            if not os.path.exists(out_file):
                include = False
                break
        if include:
            afs.append(align_file)

    print "Evaluating dataset %r: %d/%d scored alignments after minor filtering" \
            % (dataset_name, len(afs), len(align_files))
    if align_files_only:
        for af in afs:
            yield af
        return

    # Iterate through score files in dataset, per alignment.
    for align_file in afs:
        scores_cols = []
        for batchscore_id in batchscore_ids:
            sc_dir = os.path.join(ds_dir, batchscore_id)
            out_file = dataset_config.get_out_file(align_file, sc_dir)
            scores = read_batchscores(out_file)
            scores_cols.append(scores)
        alignment = Alignment(align_file,
                test_file=dataset_config.get_test_file(align_file),
                parse_testset_fn=dataset_config.parse_testset_fn)
        yield alignment, scores_cols



################################################################################
# Cmd line driver
################################################################################

def get_evaluator(name):
    try:
        ev_module = imp.load_source(name, os.path.join('evaluators', '%s.py' % name))
        fn_name = name.split('.')[-1]
        fn = getattr(ev_module, fn_name)
    except (ImportError, AttributeError), e:
        raise ImportError("%s: %s is not a valid evaluator." % (e, name))
    return fn


def main():
    parser = argparse.ArgumentParser(
        description="Run evaluators of scoring methods using the output from batch running the scorers on a dataset(s)",
        usage="%(prog)s [-h] [-l] evaluator_name dataset_name <batchscore_ids...>")

    parser.add_argument('-l', dest='list_batchscore_ids', type=str,
        help="list all batchscore run ids for the dataset")

    parser.add_argument('evaluator_name', nargs='?',
        help="name of evaluator")
    parser.add_argument('dataset_name', nargs='?',
        help="name of dataset")
    parser.add_argument('batchscore_ids', nargs='*',
        help="ids of batchscore runs to evaluate")

    parser.add_argument('-p', dest='evaluator_params', action='append', default=[],
        help="parameters to pass to the evaluator, can specify multiple. Specify as '-p paramName=paramValue', e.g. '-p overwrite=1")


    args = parser.parse_args()
    if args.list_batchscore_ids:
        ba_dir = get_batchscore_dir(args.list_batchscore_ids)
        fnames = os.listdir(ba_dir)
        print "\n".join(sorted(x for x in fnames if os.path.isdir(os.path.join(ba_dir,x))))
        sys.exit(0)
    if not args.evaluator_name or not args.dataset_name:
        parser.print_usage()
        print "%s: error: too few arguments" % sys.argv[0]
        sys.exit(1)

    args.evaluator_params = parse_params(args.evaluator_params)

    ev_name = args.evaluator_name
    ev_fn = get_evaluator(ev_name)

    ev_fn(args.dataset_name, *args.batchscore_ids, **args.evaluator_params)


if __name__ == "__main__":
    main()
