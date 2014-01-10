#!/usr/bin/python
"""
Score multiple alignments using multiple scorers.  Save scores for future
retrieval in evaluators.
"""
import argparse
import multiprocessing
import os
import sys
import time
import yaml

from conseval.alignment import Alignment
from conseval.datasets import DATASET_CONFIGS, OUTPUT_DIR
from conseval.io import write_batchscores, list_scorer_params
from conseval.scorer import get_scorer
from conseval.utils import parallelize



################################################################################
# Input/output
################################################################################

def read_batchscore_config(config_file):
    with open(config_file) as f:
        config_yaml = f.read()
    config = yaml.load(config_yaml)

    # List of dataset names
    datasets = config['datasets']

    # List of scorer id's, names, and params.
    # Initialize the scorers
    scorers = []
    scs = config['scorers']
    ids = set()
    for sc in scs:
        scorer_id = sc['id']
        if scorer_id in ids:
            raise ValueError("Duplicate scorer id: %s" % scorer_id)
        ids.add(scorer_id)
        scorer_name = sc['name']
        if 'params' in sc:
            params = sc['params']
        else:
            params = {}
        scorer = get_scorer(scorer_name, **params)
        scorer.set_output_id(scorer_id)
        scorers.append(scorer)
    return datasets, scorers



################################################################################
# Parallelization routines and helpers
################################################################################

def run_experiments(dataset_name, scorers, limit=0):
    """
    Returns iterator over lists of tuples of scores.

    The tuple is over all scoring methods in scorers.
    The list is over all sites in the alignment.
    The iterator is over all alignments in `dataset_name`.

    @param dataset_name:
        name of dataset to use
    @param scorers:
        scorers to use
    @param limit:
        Max number of alignments to score (TODO)
    """
    dataset_config = DATASET_CONFIGS[dataset_name]
    align_files = dataset_config.get_align_files(limit)

    count = 0
    tot = len(align_files)
    sys.stderr.write("\nDataset: %s\n" % dataset_name)
    sys.stderr.write("Scorers:\n")
    for scorer in scorers:
        sys.stderr.write("\t%s --> %s\n" % (scorer.output_id, scorer.output_dir))
    sys.stderr.write("%d alignments to score\n" % tot)

    if not align_files:
        return

    run_experiment = run_experiment_helper(dataset_config, scorers)

    # Shortcut if no parallelization.  Also helps debugging.
    no_parallel = (len(align_files) == 1)

    if no_parallel:
        it = ((af, run_experiment(af)) for af in align_files)
    else:
        it = parallelize.imap_unordered(run_experiment, align_files)

    t0 = time.time()
    t00 = t0
    for align_file, success in it:
        count += 1
        if time.time() - t0 > 60:
            dt = time.time() - t00
            h = dt // 3600
            m = (dt % 3600) // 60
            s = dt % 60
            sys.stderr.write("\nTime elapsed: %dh %dm %ds\n" % (h,m,s))
            sys.stderr.write("Progress: %d / %d\n" % (count, tot))
            t0 = time.time()


def run_experiment_helper(dataset_config, scorers):
    def run_experiment(align_file):
        """
        Run scorers on one aln file.  This is a helper for multithreading the
        scoring of each aln file.
        """
#       test_file = dataset_config.get_test_file(align_file)
        alignment = Alignment(align_file)
#                test_file=test_file,
#                parse_testset_fn=dataset_config.parse_testset_fn)
        for scorer in scorers:
            # Score.
            try:
                scores = scorer.score(alignment)
            except Exception, e:
                import traceback
                sys.stderr.write("\nError scoring %s via %s\n" %
                    (alignment.align_file, type(scorer).__name__))
                traceback.print_exc()
                continue
            # Write scores.
            out_file = dataset_config.get_out_file(align_file, scorer.output_dir)
            write_batchscores(out_file, scores)
        return True
    return run_experiment



################################################################################
# Cmd line driver
################################################################################

def main():
    parser = argparse.ArgumentParser(
        description="Batch score the conservation of multiple alignments, using multiple scorers.")

    parser.add_argument('config_file',
        help="YAML config file specifying the dataset and scorers.  See `exapmles/example.yaml` for an example.")
    args = parser.parse_args()


    dataset_names, scorers = read_batchscore_config(args.config_file)

    # Sanity check the output dirs
    for ds_name in dataset_names:
        ds_dir = os.path.join(OUTPUT_DIR, "batchscore-%s" % ds_name)
        if not os.path.exists(ds_dir):
            os.mkdir(ds_dir)
        for scorer in scorers:
            sc_dir = os.path.join(ds_dir, scorer.output_id)
            if os.path.exists(sc_dir):
                resp = raw_input("%s exists. Overwrite? y/[n]: " % sc_dir)
                if resp != 'y':
                    sys.exit(0)
                try:
                    for filename in os.listdir(sc_dir):
                        os.remove(os.path.join(sc_dir, filename))
                    os.rmdir(sc_dir)
                except OSError:
                    raise OSError("Could not overwrite directory %s" % sc_dir)

    # Perform the scoring
    for ds_name in dataset_names:
        ds_dir = os.path.join(OUTPUT_DIR, "batchscore-%s" % ds_name)
        for scorer in scorers:
            sc_dir = os.path.join(ds_dir, scorer.output_id)
            os.mkdir(sc_dir)
            params_file = os.path.join(ds_dir, "%s.params" % scorer.output_id)
            with open(params_file, 'w') as f:
                f.write(list_scorer_params(scorer))
            scorer.set_output_dir(sc_dir)
        run_experiments(ds_name, scorers)


if __name__ == "__main__":
    main()
