from __future__ import division
import bisect
import numpy as np
import matplotlib.pyplot as plt
import random
from sklearn.metrics import roc_curve, precision_recall_curve

from evaluate import get_batchscores, get_out_dir
from conseval.utils.stats import zscore


AUC_LEVELS = [.1, .5, 1]


_f = (.0105/.064 - .0105)
_F = 1 - .0105
def to_keep(ts):
    # need to make 1.1 pos -> 6.4 pos
    # if not pos, then keep with probability f/F
    # f = (.011/.064 - .011)
    # F = 1 - .011
    return ts or random.random() <= _f/_F


def roc_pr(dataset_name, *scorer_ids):
    """
    build roc curve for each alignment, for each scorer
    """
    scores_cols = [[] for i in scorer_ids]
    test_scores = []

    # Just aggregate all scores across all data files.  It isn't much memory anyway.
    for alignment, scores_tup in get_batchscores(dataset_name, scorer_ids):
        ts = alignment.testset
        counts = []
        for scores_col, scores in zip(scores_cols, scores_tup):
            counts.append(ts.count(None))
            scores_col += [scores[i] for i in xrange(len(ts)) if ts[i] is not None]
        test_scores += [ts[i] for i in xrange(len(ts)) if ts[i] is not None]

    print test_scores.count(1) / len(test_scores)

    scorer_fprs = []
    scorer_tprs = []
    scorer_precisions = []
    scorer_recalls = []
    for scores_col in scores_cols:
        fprs, tprs, _ = roc_curve(test_scores, scores_col, pos_label=1)
        scorer_fprs.append(fprs)
        scorer_tprs.append(tprs)
        precisions, recalls, _ = precision_recall_curve(test_scores, scores_col, pos_label=1)
        scorer_precisions.append(precisions)
        scorer_recalls.append(recalls)

    print_auc("ROC", scorer_ids, scorer_fprs, scorer_tprs)
    print_auc("PR", scorer_ids, scorer_recalls, scorer_precisions)

    fig1 = plot_roc(dataset_name, scorer_fprs, scorer_tprs, scorer_ids, .5, None)
    fig2 = plot_pr(dataset_name, scorer_precisions, scorer_recalls, scorer_ids)

#    import ipdb
#    ipdb.set_trace()

    plt.show()


def plot_roc(name, scorer_fprs, scorer_tprs, scorer_ids, x_max, y_min, legend='lower right'):
    if x_max and y_min:
        raise Exception()

    fig = plt.figure()

    x_min = 1
    y_max = 0
    for fprs, tprs, scorer_id in zip(scorer_fprs, scorer_tprs, scorer_ids):
        if x_max:
            ind = bisect.bisect(fprs, x_max)
            fprs = fprs[:ind+1]
            tprs = tprs[:ind+1]
            y_max = max(y_max, tprs[-1])
        elif y_min:
            ind = bisect.bisect(tprs, y_min)
            fprs = fprs[ind:]
            tprs = tprs[ind:]
            x_min = min(x_min, fprs[0])
        plt.plot(fprs, tprs, label=scorer_id)
    plt.plot([0,1],[0,1],'k--')

    if x_max:
        x_rg = [0,x_max]
        y_rg = [0,y_max]
    elif y_min:
        x_rg = [x_min,1]
        y_rg = [y_min,1]
    else:
        x_rg = y_rg = [0,1]
    plt.xlim(x_rg)
    plt.ylim(y_rg)
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.title('%s - ROC curve' % name)
    if legend:
        plt.legend(loc=legend)

    return fig


def plot_pr(name, scorer_precisions, scorer_recalls, scorer_ids, legend='upper right'):
    fig = plt.figure()
    y_max = 0
    for ps, rs, id in zip(scorer_precisions, scorer_recalls, scorer_ids):
        y_max = max(y_max, np.max(ps[rs>.1]))
        plt.plot(rs, ps, label=id)
    plt.xlabel('Recall')
    plt.ylabel('Precision')
    plt.ylim([0.0, y_max+.05])
    plt.xlim([0.0, 1.0])
    plt.title('%s - Precision-Recall curve' % name)
    if legend:
        plt.legend(loc=legend)
    return fig


def print_auc(name, scorer_ids, scorer_xs, scorer_ys):
    # Compute AUCs for each scorer
    scorers_aucs = []
    for xs, ys in zip(scorer_xs, scorer_ys):
        scorer_aucs = []
        for auc_level in AUC_LEVELS:
            inds_bool = xs < auc_level
            auc = np.trapz(ys[inds_bool], xs[inds_bool])
            scorer_aucs.append(auc)
        scorers_aucs.append(scorer_aucs)

    print ""
    print "%s:" % name
    print "%s\tScorer" % "\t".join("AUC_{%.2f}"%aucl for aucl in AUC_LEVELS)
    for scorer_id, scorer_auc in zip(scorer_ids, scorers_aucs):
        line = []
        for auc in scorer_auc:
            line.append("%.4f"%auc)
        line.append(scorer_id)
        print "\t".join(line)



