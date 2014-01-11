#!/usr/bin/python
import argparse
import imp
import os
import sys

from conseval.alignment import Alignment
from conseval.datasets import DATASET_CONFIGS, OUTPUT_DIR
from conseval.io import read_batchscores
from conseval.utils.bio import get_column
from conseval.utils.general import get_timestamp



def get_batchscore_dir(dataset_name):
    return os.path.join(OUTPUT_DIR, "batchscore-%s" % dataset_name)


def get_out_dir(evaluator_id=None):
    if not evaluator_id:
        evaluator_id = get_timestamp()
    ev_dir = os.path.join(OUTPUT_DIR, "evaluate-%s" % evaluator_id)
    if os.path.exists(ev_dir):
        resp = raw_input("%s exists. Overwrite? y/[n]: " % ev_dir)
        if resp != 'y':
            sys.exit(0)
        try:
            for filename in os.listdir(ev_dir):
                os.remove(os.path.join(ev_dir, filename))
            os.rmdir(ev_dir)
        except OSError:
            raise OSError("Could not overwrite directory %s" % ev_dir)
    os.mkdir(ev_dir)
    return ev_dir


def get_batchscores(dataset_name, batchscore_ids=[]):
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

    # Iterate through score files in dataset, per alignment.
    for align_file in afs:
        scores_tup = []
        for batchscore_id in batchscore_ids:
            sc_dir = os.path.join(ds_dir, batchscore_id)
            out_file = dataset_config.get_out_file(align_file, sc_dir)
            scores = read_batchscores(out_file)
            scores_tup.append(scores)
        alignment = Alignment(align_file,
                test_file=dataset_config.get_test_file(align_file),
                parse_testset_fn=dataset_config.parse_testset_fn)
        yield alignment, scores_tup



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
        description="Run evaluators of scoring methods using the output from batch running the scorers on a dataset(s)")

    parser.add_argument('evaluator_name',
        help="name of evaluator")
    parser.add_argument('dataset_name',
        help="name of dataset")
    parser.add_argument('batchscore_ids', nargs='*',
        help="ids of batchscore runs to evaluate")

    parser.add_argument('-l', dest='list_batchscore_ids', action='store_true',
        help="list all batchscore run ids for the dataset")

    args = parser.parse_args()
    if args.list_batchscore_ids:
        ds_dir = os.path.join(OUTPUT_DIR, 'batchscore-'+args.dataset_name)
        fnames = os.listdir(ds_dir)
        print "\n".join(sorted(x for x in fnames if os.path.isdir(os.path.join(ds_dir,x))))
        sys.exit(0)

    ev_name = args.evaluator_name
    ev_fn = get_evaluator(ev_name)

    ev_fn(args.dataset_name, *args.batchscore_ids)


if __name__ == "__main__":
    main()
