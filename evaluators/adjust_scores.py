import os
import sys
from conseval.datasets import DATASET_CONFIGS
from conseval.io import read_batchscores, write_batchscores
from conseval.utils.stats import norm_scores, window_scores
from evaluate import get_batchscore_dir


def adjust_scores(dataset_name, *batchscore_ids, **kwargs):
    """
    Create normalized versions of batchscore runs in `batchscore_ids`.
    Note that zscores with an absolute of higher than 5 are chopped off (set to 5).
    """
    if 'overwrite' in kwargs:
        overwrite = kwargs['overwrite']
    else:
        overwrite = False

    # Determine adjustment type
    adj_type = None
    for k,v in kwargs.iteritems():
        if k in ['norm', 'window', 'neg']:
            if adj_type:
                raise Exception("Cannot adjust in multiple ways: %s, %s" % (adj_type, k))
            if k == 'window':
                window_size = int(v)
                adj_type = '%s_%d' % (k, window_size)
            else:
                adj_type = k

    # Perform adjustments on all scores..
    for batchscore_id in batchscore_ids:
        dc = DATASET_CONFIGS[dataset_name]

        ds_dir = get_batchscore_dir(dataset_name)
        sc_dir = os.path.join(ds_dir, batchscore_id)
        sc_dir_adj = "%s-%s" % (sc_dir, adj_type)
        if not overwrite:
            resp = raw_input("Creating %s'd version of\n  %s\nat\n  %s\nContinue? y/[n]: "
                    % (adj_type, sc_dir, sc_dir_adj))
            if resp != 'y':
                continue

        if not os.path.exists(sc_dir_adj):
            os.mkdir(sc_dir_adj)
        else:
            for fname in os.listdir(sc_dir_adj):
                os.remove(os.path.join(sc_dir_adj, fname))

        align_files = dc.get_align_files()
        for align_file in align_files:
            out_file = dc.get_out_file(align_file, sc_dir)
            if not os.path.exists(out_file):
                continue
            scores = read_batchscores(out_file)

            if adj_type is 'norm':
                scores_adj = norm_scores(scores, filter=5)
            elif adj_type is 'neg':
                scores_adj = [-s for s in scores]
            elif adj_type.startswith('window'):
                scores_adj = window_scores(scores, window_size)

            out_file_adj = os.path.join(sc_dir_adj, os.path.split(out_file)[-1])
            write_batchscores(out_file_adj, scores_adj)

        print "Created %s" % sc_dir_adj
