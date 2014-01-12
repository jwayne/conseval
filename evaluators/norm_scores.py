import os
import sys
from conseval.datasets import DATASET_CONFIGS
from conseval.io import read_batchscores, write_batchscores
from conseval.utils.stats import zscore
from evaluate import get_batchscore_dir


def norm_scores(dataset_name, *batchscore_ids, **kwargs):
    """
    Create normalized versions of batchscore runs in `batchscore_ids`.
    Note that zscores with an absolute of higher than 5 are chopped off (set to 5).
    """
    if 'overwrite' in kwargs:
        overwrite = kwargs['overwrite']
    else:
        overwrite = False

    for batchscore_id in batchscore_ids:
        dc = DATASET_CONFIGS[dataset_name]

        ds_dir = get_batchscore_dir(dataset_name)
        sc_dir = os.path.join(ds_dir, batchscore_id)
        sc_dir_norm = sc_dir + "-norm"
        if not overwrite:
            resp = raw_input("Creating normalized version of\n  %s\nat\n  %s\nContinue? y/[n]: "
                    % (sc_dir, sc_dir_norm))
            if resp != 'y':
                continue

        if not os.path.exists(sc_dir_norm):
            os.mkdir(sc_dir_norm)
        else:
            for fname in os.listdir(sc_dir_norm):
                os.remove(os.path.join(sc_dir_norm, fname))

        align_files = dc.get_align_files()
        for align_file in align_files:
            out_file = dc.get_out_file(align_file, sc_dir)
            if not os.path.exists(out_file):
                continue
            scores = read_batchscores(out_file)

            scores_norm = zscore(scores, filter=5)
            out_file_norm = os.path.join(sc_dir_norm, os.path.split(out_file)[-1])
            write_batchscores(out_file_norm, scores_norm)

        print "Created %s" % sc_dir_norm
