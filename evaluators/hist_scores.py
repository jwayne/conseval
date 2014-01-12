from __future__ import division
import matplotlib.pyplot as plt
import numpy as np
from evaluate import get_batchscores


def hist_scores(dataset_name, *batchscore_ids, **kwargs):
    """
    Histogram positive and negative scores for each scorer.  Plot F1 for corresponding
    thresholds alongside.
    """
    if 'noblock' in kwargs:
        noblock = kwargs['noblock']
    else:
        noblock = False

    N = len(batchscore_ids)
    pos_cols = [[] for i in xrange(N)]
    neg_cols = [[] for i in xrange(N)]

    # Just aggregate all scores across all data files.  It isn't much memory anyway.
    for alignment, scores_cols in get_batchscores(dataset_name, batchscore_ids):
        ts = alignment.testset
        for pos_col, neg_col, scores_col in zip(pos_cols, neg_cols, scores_cols):
            pos_col += (scores_col[i] for i in xrange(len(scores_col)) if ts[i])
            neg_col += (scores_col[i] for i in xrange(len(scores_col)) if not ts[i])

    figs = []
    for pos_col, neg_col, batchscore_id in zip(pos_cols, neg_cols, batchscore_ids):
        pos_col = np.array(pos_col)
        neg_col = np.array(neg_col)
        scores_min = np.min((np.min(pos_col), np.min(neg_col)))
        scores_max = np.max((np.max(pos_col), np.max(neg_col)))
        bins = np.linspace(scores_min,scores_max,101)
        bin_width = bins[1] - bins[0]
        bin_lows = bins[:-1]
        bin_mids = bin_lows + bin_width

        # Plot counts
        fig = plt.figure()
        pos_counts, _ = np.histogram(pos_col, bins)
        neg_counts, _ = np.histogram(neg_col, bins)
        ax = plt.gca()
        ax.bar(bin_lows, pos_counts, bin_width, color='g', alpha=0.8)
        ax.bar(bin_lows, neg_counts, bin_width, color='r', bottom=pos_counts, alpha=0.5)
        ax.set_ylabel('Count')

        # Plot corresponding F1
        pos_tot = sum(pos_counts)
        pos_right = 0
        neg_right = 0
        f1 = []
        for pos_count, neg_count in reversed(zip(pos_counts, neg_counts)):
            pos_right += pos_count
            neg_right += neg_count
            precision = pos_right / pos_tot
            recall = pos_right / (pos_right + neg_right)
            if precision or recall:
                f1.append( 2 * precision * recall / (precision + recall) )
            else:
                f1.append(0)
        f1 = np.array(list(reversed(f1)))
        ax = ax.twinx()
        ax.plot(bin_mids, f1, '--')
        ax.set_ylabel('F1')

        plt.title('Scores (green=+, red=-, --=F1) on %s by %s' % (dataset_name, batchscore_id))
        plt.xlabel('Score')
        figs.append(fig)

    if noblock:
        plt.show()
    else:
        plt.show(block=False)
        import ipdb
        ipdb.set_trace()
