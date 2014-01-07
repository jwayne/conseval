import os
import sys
from conseval.datasets import DATASET_CONFIGS
from conseval.io import read_batchscores, write_batchscores
from evaluate import get_batchscore_dir


def convert_scores_neg(dataset_name, *scorer_ids):
    for scorer_id in scorer_ids:
        dc = DATASET_CONFIGS[dataset_name]

        ds_dir = get_batchscore_dir(dataset_name)
        sc_dir = os.path.join(ds_dir, scorer_id)
        sc_dir_neg = sc_dir + "-neg"
        resp = raw_input("Creating negated version of\n  %s\nat\n  %s\nContinue? y/[n]: "
                % (sc_dir, sc_dir_neg))
        if resp != 'y':
            continue

        if not os.path.exists(sc_dir_neg):
            os.mkdir(sc_dir_neg)

        align_files = dc.get_align_files()
        for align_file in align_files:
            out_file = dc.get_out_file(align_file, sc_dir)
            scores = read_batchscores(out_file)

            scores_neg = [-s for s in scores]
            out_file_neg = os.path.join(sc_dir_neg, os.path.split(out_file)[-1])
            write_batchscores(out_file_neg, scores_neg)
