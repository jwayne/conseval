#! /usr/bin/python

import numpy as np
import os
import sys

from alignment import Alignment
from process import run_scorers
from scorer import get_scorer
from utils import get_column


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

def run_experiments(dataset='all', scorer_names='js_divergence', **kwargs):
    # Determine which sets of data to run experiments on
    if dataset == 'all':
        data_configs = DATA_CONFIGS.values()
    elif type(dataset) == str:
        data_configs = (DATA_CONFIGS[dataset],)
    else:
        data_configs = [DATA_CONFIGS[ds] for ds in dataset]

    # Determine which scorer to run
    if type(scorer_names) == str:
        scorer_names = (scorer_names,)
    scorers = [get_scorer(s) for s in scorer_names]

    count_scored = 0
    count_missing = 0
    for data_config in data_configs:
        aln_dir = os.path.join(DATA_HOME_DIR, data_config['aln_dir'])
        for root, dirs, files in os.walk(aln_dir):
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
                count_scored += 1

                alignment = Alignment(align_file)
                scores = run_scorers(alignment, scorers, **kwargs)
                actual = parse_testset(test_file, data_config['parse_fields_func'], alignment)

                yield [(x,y) for x,y in zip(scores, actual) if y is not None]
    sys.stderr.write("Scored %d alignments, missing testfiles for %d alignments\n"
        % (count_scored, count_missing))


if __name__ == "__main__":
    for exp in run_experiments():
        print exp[0]
        break
