#!/usr/bin/python
import argparse
import multiprocessing
import os
import random
import shutil
import sys
import yaml

from alignment import Alignment
from datasets import DATASET_CONFIGS, OUTPUT_DIR
from scorers.base import get_scorer
from singlerun import compute_scores, prepare_header, write_scores, read_scores
from utils import parallelize
from utils.general import get_timestamp


def run_experiments(scorers, dataset_name, out_dir, limit=0, section=None):
    """
    Returns iterator over lists of tuples of scores.

    The tuple is over all scoring methods in scorers.
    The list is over all sites in the alignment.
    The iterator is over all alignments in `dataset_name`.

    @param scorers:
        scorers to use
    @param dataset_name:
        name of dataset to use
    @param out_dir:
        write output of scorers to this directory
    @param limit:
        Max number of alignments to score
    @param section:
        section of the dataset to score
    """
    # Handle dataset_name, limit
    dataset_config = DATASET_CONFIGS[dataset_name]
    align_files = dataset_config.get_align_files(limit, section)
    tot = len(align_files)
    sys.stderr.write("Scoring %d alignments to %s\n" % (tot, out_dir))
    if not align_files:
        return

    run_experiment = run_experiment_helper(scorers, dataset_config, out_dir)

    # Shortcut if no parallelization.  Also helps debugging.
    no_parallel = (len(align_files) == 1)

    if no_parallel:
        it = ((af, run_experiment(af)) for af in align_files)
    else:
        it = parallelize.imap_unordered(run_experiment, align_files)

    for align_file, score_tups in it:
        pass


def run_experiment_helper(scorers, dataset_config, out_dir):
    scorer_names = [scorer.name for scorer in scorers]
    def run_experiment(align_file):
        """
        Run scorers on one aln file.  This is a helper for multithreading the
        scoring of each aln file.
        """
#       test_file = dataset_config.get_test_file(align_file)
        alignment = Alignment(align_file)
#                test_file=test_file,
#                parse_testset_fn=dataset_config.parse_testset_fn)
        score_tups = compute_scores(alignment, scorers)
        out_file = dataset_config.get_out_file(align_file, out_dir)
        with open(out_file, 'w') as f:
            # TODO: we don't really need to output the first 2 columns in
            # write_scores when in batch mode, since usually when we read the
            # files for evaluations we ignore the first 2 columns (i.e. pos and site)
            write_scores(alignment, score_tups, scorer_names, f)
        return score_tups
    return run_experiment


################################################################################
# Input/output
################################################################################

def read_config_file(config_file):
    with open('config.yaml') as f:
        config_yaml = f.read()
    config = yaml.load(config_yaml)
    return config

def read_config_datasets(config):
    # Collect dataset names/parameters in config file
    datasets = []
    das = config['datasets']
    ids = set()
    for da in das:
        dataset_id = da['id']
        if dataset_id in ids:
            raise ValueError("Duplicate dataset id: %s" % dataset_id)
        ids.add(dataset_id)
        dataset_name = da['name']
        dc_params = {}
        if 'limit' in da:
            dc_params['limit'] = da['limit']
        if 'section' in da:
            dc_params['section'] = da['section']
        datasets.append((dataset_id, dataset_name, dc_params))
    return datasets

def read_config_scorers(config):
    # Initialize scorers in config file
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
        scorer.name = scorer_id
        scorers.append(scorer)
    return scorers


################################################################################
# Cmd line driver
################################################################################

def main():
    parser = argparse.ArgumentParser(
        description="Produce conservation scores for a dataset via multiple scorers.")

    parser.add_argument('config_file',
        help="YAML config file specifying the dataset and scorers")
    args = parser.parse_args()


    config = read_config_file(args.config_file)
    datasets = read_config_datasets(config)
    scorers = read_config_scorers(config)

    # Make the output directory.
    count = 0
    batch_dir = os.path.join(OUTPUT_DIR, "batch-%s-0" % config['id'])
    while os.path.exists(batch_dir):
        count += 1
        batch_dir = os.path.join(OUTPUT_DIR, "batch-%s-%d" % (config['id'], count))
    os.mkdir(batch_dir)
    sys.stderr.write("Output directory: %s\n" % batch_dir)

    # Create param/config files in output directory.
    with open(os.path.join(batch_dir, "params.txt"), 'w') as f:
        f.write(prepare_header(scorers))
    shutil.copy(args.config_file, os.path.join(batch_dir, "config.yaml"))

    # Run scorers on each dataset, outputting to separate dirs for different datasets
    for dataset_id, dataset_name, dc_params in datasets:
        out_dir = os.path.join(batch_dir, dataset_id)
        os.mkdir(out_dir)
        run_experiments(scorers, dataset_name, out_dir, **dc_params)


if __name__ == "__main__":
    main()
