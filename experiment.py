#!/usr/bin/python
import argparse
import datetime
import multiprocessing
import numpy as np
import os
import random
import sys

from alignment import Alignment
from scorer import get_scorer
from singlerun import compute_scores, write_scores, parse_scorer_names, DEFAULT_SCORER
from utils.bio import get_column
from utils import parallelize


################################################################################
# Data processing
################################################################################

def parse_testset(test_file, parse_fields_func, alignment):
    actual = []
    start_pos = None

    with open(test_file, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            fields = line.split()
            ind, result = parse_fields_func(fields)
            if start_pos is None:
                start_pos = ind
            pos = ind - start_pos

            if len(actual) > pos:
                # if indices skip down, then re-adjust so it's normal
                start_pos -= len(actual) - pos
            elif pos > len(actual):
                for i in xrange(len(actual), pos):
                    if get_column(i, alignment.msa) == 'X':
                        # if shitty, fill it in
                        actual.append(None)
                        sys.stderr.write("%d: %s\n" % (ind, test_file))
                    else:
                        # otherwise, re-adjust
                        start_pos += 1
            actual.append(result)
    return actual

def parse_testset_csa_fields(fields):
    #result is 1 if catalytic site, 0 if not
    result = int(int(fields[1]) > 0)
    pos = int(fields[0])
    return pos, result

def parse_testset_ec_fields(fields):
    #result is distance of site from ligand atom
    result = float(fields[2])
    pos = int(fields[0])
    return pos, result


DATA_HOME_DIR = os.path.abspath("../../data/cs08/")

DATA_CONFIGS = {
    'csa': {
        'aln_dir': 'conservation_alignments/csa_hssp',
        'test_dir': 'cat_sites',
        'aln_to_test': lambda x: x.split('.')[0] + '.cat_sites',
        'parse_fields_func': parse_testset_csa_fields,
    },
    'ec': {
        'aln_dir': 'conservation_alignments/ec_hssp',
        'test_dir': 'lig_distance',
        'aln_to_test': lambda x: x[:6]+ '.dist_to_lig',
        'parse_fields_func': parse_testset_ec_fields,
    },
}


################################################################################
# Run experiments
################################################################################

def run_experiments(data_names='all', scorer_config={'js_divergence':{}},
        limit=0, out_dirname=None):
    """
    Iterator that returns lists of pairs of observed/expected scores for each
    column in a file, for all alignment files in the datasets requested.  Each
    item in the iterator has the form:
        [( (score1, score2, score3, ...), expected_score ) for col1,
         ( (score1, score2, score3, ...), expected_score ) for col2,
         ...
        ]
    @param data_names:
        'all', dataset to use, or list of datasets to use
    @param scorer_config:
        Dict mapping scorer names to their desired parameters
    @param limit:
        Max number of alignments to score
    """
    # Determine which sets of data to run experiments on
    if data_names == 'all':
        data_configs = DATA_CONFIGS.values()
    elif type(data_names) == str:
        data_configs = (DATA_CONFIGS[data_names],)
    else:
        data_configs = [DATA_CONFIGS[ds] for ds in data_names]

    # Determine which scorer to run
    scorers = []
    for scorer_name in sorted(scorer_config.keys()):
        scorer_params = scorer_config[scorer_name]
        scorers.append(get_scorer(scorer_name, **scorer_params))

    if out_dirname:
        sys.stderr.write("Writing scores to %s/\n" % os.path.abspath(out_dirname))

    # Assemble list of files to score.
    filenames = get_filenames(data_configs, limit)

    run_experiment = run_experiment_helper(scorers, out_dirname)

    # Shortcut if no parallelization
    if len(filenames) == 1:
        yield run_experiment(filenames[0])
        return

    it = parallelize.imap_unordered(run_experiment, filenames)
    for arg, result in it:
        align_file = arg[0]
        all_scores = result
        yield align_file, all_scores


def get_filenames(data_configs, limit):
    count_scoring = 0
    count_missing = 0
    filenames = []
    for data_config in data_configs:
        aln_dir = os.path.join(DATA_HOME_DIR, data_config['aln_dir'])
        for root, _, files in os.walk(aln_dir):
            for file in files:
                if not file.endswith('.aln'):
                    continue
                align_file = os.path.join(root, file)
                # Get the 'tail' of the path after aln_dir
                file = os.path.join(root, data_config['aln_to_test'](file))[len(aln_dir)+1:]
                test_file = os.path.join(DATA_HOME_DIR, data_config['test_dir'], file)
                if not os.path.exists(test_file):
                    # No test file, so don't even bother scoring.
                    count_missing += 1
                    continue
                filenames.append((align_file, test_file, data_config['parse_fields_func']))
                count_scoring += 1
    sys.stderr.write("Found %d alignments w/ testfiles, %d alignments w/o testfiles\n"
            % (count_scoring, count_missing))
    if limit and len(filenames) > limit:
        random.seed(1000)
        filenames = random.sample(filenames, limit)
    sys.stderr.write("Scoring %d alignments\n" % len(filenames))
    return filenames


def run_experiment_helper(scorers, out_dir):
    scorer_names = [scorer.name for scorer in scorers]
    def run_experiment(args):
        """
        Run scorers on one aln file.  This is a helper for multithreading the
        scoring of each aln file.
        """
        align_file, test_file, parse_fields_func = args
        alignment = Alignment(align_file)
        all_scores = [parse_testset(test_file, parse_fields_func, alignment)]
        score_tups = compute_scores(alignment, scorers, all_scores)
        if out_dir:
            out_fname = ".".join(os.path.split(align_file)[-1].split('.')[:-1]) + ".res"
            with open(os.path.join(out_dir, out_fname), 'w') as f:
                write_scores(alignment, score_tups, scorer_names, f,
                        includes_testset=True)
    return run_experiment




################################################################################
# Cmd line driver
################################################################################

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run experiments of computing conservation scores on a set of alignment files and comparing those results with the test data.")
    parser.add_argument('out_dirname')

    # TODO: these settings suck
    parser.add_argument('-d', dest='data_names', action='append', default=[],
        help="datasets to run experiments on.")
    parser.add_argument('-s', dest='scorer_names', action='append', default=[DEFAULT_SCORER],
        help="conservation estimation methods")
    parser.add_argument('-n', dest='limit', type=int, default=0,
        help="max number of alignment files to run the experiments on.")
    args = parser.parse_args()
    if not args.data_names:
        args.data_names = 'all'
    elif len(args.data_names) == 1:
        args.data_names = args.data_names[0]
    scorer_config = dict((name, {}) for name in args.scorer_names)

    ts = datetime.datetime.now().strftime("%Y%m%d-%H%M%S")
    out_dirname = os.path.join(args.out_dirname, "experiment-%s" % ts)
    os.mkdir(out_dirname)

    try:
        for align_file, all_scores in run_experiments(
                data_names=args.data_names,
                scorer_config=scorer_config,
                limit=args.limit,
                out_dirname=out_dirname):
            pass
    finally:
        if not os.listdir(out_dirname):
            os.rmdir(out_dirname)
