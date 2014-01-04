#!/usr/bin/python
import argparse
import datetime
import multiprocessing
import numpy as np
import os
import random
import sys
import yaml

from alignment import Alignment
from datasets import DATASET_CONFIGS, OUT_HOME_DIR
from scorers.base import get_scorer
from single import compute_scores, prepare_header, write_scores
from utils.bio import get_column
from utils import parallelize
from utils.general import get_timestamp


################################################################################
# Run experiments
################################################################################

def run_experiments(scorer_config, dataset_name, limit=0, no_parallel=False):
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
    ts = get_timestamp()
    out_dirname = os.path.join(OUT_HOME_DIR, "batch-%s" % ts)
    sys.stderr.write("Writing scores to %s/\n" % out_dirname)
    os.mkdir(out_dirname)

    # Print header of params used in this experiment.
    with open(os.path.join(out_dirname, "_params.txt"), 'w') as f:
        f.write(prepare_header(scorers))

    run_experiment = run_experiment_helper(scorers, dataset_config, out_dirname)

    # Shortcut if no parallelization.  Also helps debugging.
    if len(align_files) == 1:
        no_parallel = True

    if no_parallel:
        it = ((af, run_experiment(af)) for af in align_files)
    else:
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
        test_file = dataset_config.get_test_file(align_file)
        alignment = Alignment(align_file, test_file=test_file,
                parse_testset_fn=dataset_config.parse_testset_fn)
        score_tups = compute_scores(alignment, scorers)
        out_file = dataset_config.get_out_file(align_file, out_dir)
        with open(out_file, 'w') as f:
            write_scores(alignment, score_tups, scorer_names, f)
        return score_tups
    return run_experiment




################################################################################
# Cmd line driver
################################################################################

def main():
    parser = argparse.ArgumentParser(
        description="Produce conservation scores for a dataset via multiple scorers.")

    parser.add_argument('config_file',
        help="YAML config file specifying the dataset and scorers")

#    parser.add_argument('-n', dest='limit', type=int, default=0,
#        help="max number of alignment files to run the experiments on.")
#    parser.add_argument('--no_parallel', action='store_true',
#        help="do not parallelize the experiments. Useful for debugging.")
    args = parser.parse_args()

    with open(args.config_file) as f:
        config_yaml = f.read()
    config = yaml.load(config_yaml)
    config[1]

    for align_file, score_tups in run_experiments(
            scorer_config=scorer_config,
            dataset_name=args.dataset_name,
            limit=args.limit):
        pass


if __name__ == "__main__":
    main()
