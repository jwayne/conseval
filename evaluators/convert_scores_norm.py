import os
import sys
from conseval.datasets import DATASET_CONFIGS
from conseval.io import read_batchscores, write_batchscores
from conseval.utils.stats import zscore
from evaluate import get_batchscore_dir


def convert_scores_norm(dataset_name, *scorer_ids):
    for scorer_id in scorer_ids:
        dc = DATASET_CONFIGS[dataset_name]

        ds_dir = get_batchscore_dir(dataset_name)
        sc_dir = os.path.join(ds_dir, scorer_id)
        sc_dir_norm = sc_dir + "-norm"
        resp = raw_input("Creating normalized version of\n  %s\nat\n  %s\nContinue? y/[n]: "
                % (sc_dir, sc_dir_norm))
        if resp != 'y':
            continue

        if not os.path.exists(sc_dir_norm):
            os.mkdir(sc_dir_norm)

        align_files = dc.get_align_files()
        for align_file in align_files:
            out_file = dc.get_out_file(align_file, sc_dir)
            if not os.path.exists(out_file):
                continue
            scores = read_batchscores(out_file)

            scores_norm = zscore(scores)
            out_file_norm = os.path.join(sc_dir_norm, os.path.split(out_file)[-1])
            write_batchscores(out_file_norm, scores_norm)
