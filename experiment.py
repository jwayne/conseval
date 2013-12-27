#!/usr/bin/python
import argparse
import datetime
import multiprocessing
import numpy as np
import os
import random
import sys

from alignment import Alignment
from dataset_config import DATASET_CONFIGS
from scorer import get_scorer
from singlerun import compute_scores, prepare_header, write_scores
from utils.bio import get_column
from utils import parallelize
from utils.general import get_timestamp


################################################################################
# Run experiments
################################################################################

def run_experiments(scorer_config, dataset_name, out_dirname, limit=0):
    """
    Iterator that returns lists of pairs of observed/expected scores for each
    column in a file, for all alignment files in the datasets requested.  Each
    item in the iterator has the form:
        [( (score1, score2, score3, ...), expected_score ) for col1,
         ( (score1, score2, score3, ...), expected_score ) for col2,
         ...
        ]
    @param scorer_config:
        Dict mapping scorer names to their desired parameters
    @param dataset_name:
        name of dataset to use
    @param out_dirname:
        directory to write scores to
    @param limit:
        Max number of alignments to score
    """
    # Handle dataset_name, limit
    dataset_config = DATASET_CONFIGS[dataset_name]
    align_files = dataset_config.get_align_files(limit)
    sys.stderr.write("Scoring %d alignments\n" % len(align_files))
    if not align_files:
        return

    # Handle scorer_config
    scorers = []
    for scorer_name in sorted(scorer_config.keys()):
        scorer_params = scorer_config[scorer_name]
        scorers.append(get_scorer(scorer_name, **scorer_params))

    # Handle out_dirname
    out_dirname = os.path.abspath(out_dirname)
    if not os.path.exists(out_dirname):
        raise IOError("Directory '%s' does not exist." % out_dirname)
    ts = get_timestamp()
    out_dirname = os.path.join(args.out_dirname, "experiment-%s" % ts)
    sys.stderr.write("Writing scores to %s/\n" % out_dirname)
    os.mkdir(out_dirname)

    # Print header of params used in this experiment.
    with open(os.path.join(out_dirname, "params.txt"), 'w') as f:
        f.write(prepare_header(scorers))

    run_experiment = run_experiment_helper(scorers, dataset_config, out_dirname)

    # Shortcut if no parallelization
    if len(align_files) == 1:
        yield align_files[0], run_experiment(align_files[0])
        return

    it = parallelize.imap_unordered(run_experiment, align_files)
    for align_file, score_tups in it:
        yield align_file, score_tups


def run_experiment_helper(scorers, dataset_config, out_dir):
    scorer_names = [scorer.name for scorer in scorers]
    def run_experiment(align_file):
        """
        Run scorers on one aln file.  This is a helper for multithreading the
        scoring of each aln file.
        """
        test_file = dataset_config.get_test_file(align_fie)
        alignment = Alignment(align_file, test_file=test_file,
                parse_testset_line=dataset_config.parse_testset_line)
#TODO: build ppc (probably just drawing distribution of rates for conserved/unconserved) for mayrose04
        score_tups = compute_scores(alignment, scorers)
        out_file = ".".join(align_file[len(dataset_config.aln_dir+1):].replace('/', '___').split('.')[:-1]) + ".res"
        with open(os.path.join(out_dir, out_file), 'w') as f:
            write_scores(alignment, score_tups, scorer_names, f)
        return score_tups
    return run_experiment




################################################################################
# Cmd line driver
################################################################################

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run experiments of computing conservation scores on a set of alignment files and comparing those results with the test data.")

    parser.add_argument('scorer_names')
    parser.add_argument('dataset_name')
    parser.add_argument('out_dirname')

    parser.add_argument('-n', dest='limit', type=int, default=0,
        help="max number of alignment files to run the experiments on.")
    args = parser.parse_args()

    scorer_names = args.scorer_names.split(',')
    scorer_config = dict((name, {}) for name in scorer_names)

    for align_file, score_tups in run_experiments(
            scorer_config=scorer_config,
            dataset_name=args.dataset_name,
            out_dirname=args.out_dirname,
            limit=args.limit):
        pass
