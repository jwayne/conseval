from __future__ import division
import matplotlib.pyplot as plt
import numpy as np
import scipy.stats as ss
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
    if 'fit_gamma' in kwargs:
        fit_gamma = kwargs['fit_gamma']
    else:
        fit_gamma = False

    batchscore_ids = list(batchscore_ids)
    N = len(batchscore_ids)
    pos_cols = [[] for i in xrange(N)]
    neg_cols = [[] for i in xrange(N)]

    # Just aggregate all scores across all data files.  It isn't much memory anyway.
    for alignment, scores_cols in get_batchscores(dataset_name, batchscore_ids):
        ts = alignment.testset
        for pos_col, neg_col, scores_col in zip(pos_cols, neg_cols, scores_cols):
            pos_col += (scores_col[i] for i in xrange(len(scores_col)) if ts[i])
            neg_col += (scores_col[i] for i in xrange(len(scores_col)) if not ts[i])

    for i in xrange(len(pos_cols)):
        pos_col = pos_cols[i]
        if isinstance(pos_col, tuple) and len(pos_col[0]) == 2:
            # We are analyzing r4s_func.
            pos_cols[i], pos_col_new = zip(*pos_col)
            pos_cols.append(pos_col_new)
            neg_cols[i], neg_col_new = zip(*neg_cols[i])
            neg_cols.append(neg_col_new)
            batchscore_ids.append(batchscore_ids[i] + '-c')

    figs = []
    for pos_col, neg_col, batchscore_id in zip(pos_cols, neg_cols, batchscore_ids):
        pos_col = np.array(pos_col)
        neg_col = np.array(neg_col)

        # Compute summary statistics
        tot = len(pos_col) + len(neg_col)
        pos_mean = np.mean(pos_col)
        pos_var = np.var(pos_col)
        neg_mean = np.mean(neg_col)
        neg_var = np.var(neg_col)
        print "%s:\n\tpos %f, neg %s\n\tpos mean %f, var %f\n\tneg mean %f, var %f" % (
            batchscore_id, len(pos_col)/tot, len(neg_col)/tot,
            pos_mean, pos_var, neg_mean, neg_var)

        # Compute plot statistics
        scores_min = np.min((np.min(pos_col), np.min(neg_col)))
        scores_max = np.max((np.max(pos_col), np.max(neg_col)))
        bins = np.linspace(scores_min,scores_max,101)
        bin_width = bins[1] - bins[0]
        bin_lows = bins[:-1]
        bin_mids = bin_lows + bin_width

        # Plot counts
        plt.ion()
        fig = plt.figure()
        plt.xlabel('Score')
        pos_counts, _ = np.histogram(pos_col, bins)
        neg_counts, _ = np.histogram(neg_col, bins)
        ax = plt.gca()
        ax.bar(bin_lows, neg_counts, bin_width, color='r', alpha=0.5)#, bottom=pos_counts)
        ax.bar(bin_lows, pos_counts, bin_width, color='g', alpha=0.8)
        ax.set_ylabel('Count')

        # Plot gamma fit if desired
        if fit_gamma:
            x = np.abs(bin_mids)
            plt.plot(bin_mids, bin_width * len(pos_col) * \
                ss.gamma.pdf(x, a=abs(pos_mean)**2/pos_var, scale=pos_var/abs(pos_mean)),
                'g')
            plt.plot(bin_mids, bin_width * len(neg_col) * \
                ss.gamma.pdf(x, a=abs(neg_mean)**2/neg_var, scale=neg_var/abs(neg_mean)),
                'r')

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

        plt.title('Scores (%s, %s)' % (batchscore_id, dataset_name))
        figs.append(fig)

    if noblock:
        plt.show()
    else:
        plt.show(block=False)
        import ipdb
        ipdb.set_trace()
