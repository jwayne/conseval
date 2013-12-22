import argparse
import multiprocessing
import numpy as np
import os
import random
import sys

from alignment import Alignment
from scorer import get_scorer, parse_scorer_names
from utils.bio import get_column


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

def run_experiments(data_config='all', scorer_config={'js_divergence':{}}, limit=0):
    """
    Iterator that returns lists of pairs of observed/expected scores for each
    column in a file, for all alignment files in the datasets requested.  Each
    item in the iterator has the form:
        [( (score1, score2, score3, ...), expected_score ) for col1,
         ( (score1, score2, score3, ...), expected_score ) for col2,
         ...
        ]
    @param data_config:
        'all', dataset to use, or list of datasets to use
    @param scorer_config:
        Dict mapping scorer names to their desired parameters
    @param limit:
        Max number of alignments to score
    """
    # Determine which sets of data to run experiments on
    if data_config == 'all':
        data_configs = DATA_CONFIGS.values()
    elif type(data_config) == str:
        data_configs = (DATA_CONFIGS[data_config],)
    else:
        data_configs = [DATA_CONFIGS[ds] for ds in data_config]

    # Determine which scorer to run
    scorers = []
    for scorer_name in sorted(scorer_config.keys()):
        scorer_params = scorer_config[scorer_name]
        scorers.append(get_scorer(scorer_name, **scorer_params))

    # Assemble list of files to score.
    filenames = get_filenames(data_configs, limit)
    args = (fn + (scorers,) for fn in filenames)

    # Shortcut if no parallelization
    if len(filenames) == 1:
        yield run_experiment(args.next())
        return

    it = parallelize(run_experiment, args)
    for result in it:
        yield result


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
    if len(filenames) > limit:
        random.seed(1000)
        filenames = random.sample(filenames, limit)
    sys.stderr.write("Scoring %d alignments" % len(filenames))
    return filenames


def run_experiment(args):
    """
    Run scorers on one aln file.  This is a helper for multithreading the
    scoring of each aln file.
    """
    align_file, test_file, parse_fields_func, scorers = args
    alignment = Alignment(align_file)
    all_scores = []
    for scorer in scorers:
        try:
            scores = scorer.score(alignment)
        except Exception, e:
            import traceback
            sys.stderr.write("Error scoring %s via %s\n" %
                (alignment.align_file, type(scorer).__name__))
            traceback.print_exc()
            continue
        all_scores.append(scores)
    actual = parse_testset(test_file, parse_fields_func, alignment)
    return [(x,y) for x,y in zip(all_scores, actual) if y is not None]




################################################################################
# For testing only
################################################################################

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run experiments of computing conservation scores on a set of alignment files and comparing those results with the test data.")
    parser.add_argument('-n', dest='limit', type=int, default=1,
        help="max number of alignment files to run the experiments on.")
    parser.add_argument('-s', dest='scorer_names', action='append', default=[],
        help="conservation estimation method")

    args = parser.parse_args()
    scorer_names = parse_scorer_names(args.scorer_names)

    for exp in run_experiments(scorer_names=scorer_names, limit=args.limit):
        import ipdb
        ipdb.set_trace()
        print exp
