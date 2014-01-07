#!/usr/bin/python
import argparse
import os
import sys

from alignment import Alignment
from batchrun import read_config_file, read_config_datasets
from datasets import DATASET_CONFIGS
from singlerun import read_scores


def get_evaluator(name):
    try:
        ev_module = __import__('evaluators.'+name, fromlist=['x'])
        fn_name = name.split('.')[-1]
        fn = getattr(ev_module, fn_name)
    except (ImportError, AttributeError), e:
        raise ImportError("%s: %s is not a valid evaluator." % (e, name))
    return fn


def _align_file_iterator(dataset_config, align_files, dataset_dir):
    for align_file in align_files:
        test_file = dataset_config.get_test_file(align_file)
        out_file = dataset_config.get_out_file(align_file, dataset_dir)
        alignment = Alignment(align_file, test_file=test_file,
                parse_testset_fn=dataset_config.parse_testset_fn)
        score_tups = read_scores(out_file)
        yield (alignment, score_tups)   



################################################################################
# Input/output
################################################################################

def read_evaluate_config(config_file):
    with open('config.yaml') as f:
        config_yaml = f.read()
    config = yaml.load(config_yaml)

    # TODO start
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
    # TODO end
    return datasets



################################################################################
# Cmd line driver
################################################################################

def main():
    parser = argparse.ArgumentParser(
        description="Run evaluators of scoring methods using the output from batch running the scorers on a dataset(s)")

    parser.add_argument('evaluator_name',
        help="name of evaluator")
    parser.add_argument('batch_output_dir',
        help="directory containing scores, i.e. output directory from batch running scorers")

    args = parser.parse_args()

    ev_name = args.evaluator_name
    ev_fn = get_evaluator(ev_name)
    batch_dir = os.path.abspath(args.batch_output_dir)

    config_file = os.path.join(batch_dir, 'config.yaml')
    config = read_config_file(config_file)
    datasets = read_config_datasets(config)
    scorer_names = [sc['id'] for sc in config['scorers']]

    for dataset_id, dataset_name, dc_params in datasets:
        dataset_dir = os.path.join(batch_dir, dataset_id)
        dc = DATASET_CONFIGS[dataset_name]
        align_files = dc.get_align_files(**dc_params)
        it = _align_file_iterator(dc, align_files, dataset_dir)
        ev_fn(it, scorer_names, batch_dir)


if __name__ == "__main__":
    main()
