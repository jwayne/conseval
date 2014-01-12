from __future__ import division
import bisect
import numpy as np
import matplotlib.pyplot as plt
import random
from sklearn.metrics import roc_curve, precision_recall_curve

from evaluate import get_batchscores
from conseval.utils.stats import zscore


# This code is temporary
_f = (.0105/.064 - .0105)
_F = 1 - .0105
def to_keep(ts):
    # need to make 1.1 pos -> 6.4 pos
    # if not pos, then keep with probability f/F
    # f = (.011/.064 - .011)
    # F = 1 - .011
    return ts or random.random() <= _f/_F
# End temporary code


def pr_roc(dataset_name, *batchscore_ids):
    """
    Draw PR and ROC curves for each scorer.
    """
    allscores_cols = [[] for i in batchscore_ids]
    test_scores = []

    # Just aggregate all scores across all data files.  It isn't much memory anyway.
    for alignment, scores_cols in get_batchscores(dataset_name, batchscore_ids):
        ts = alignment.testset
        counts = []
        for allscores_col, scores in zip(allscores_cols, scores_cols):
            counts.append(ts.count(None))
            allscores_col += [scores[i] for i in xrange(len(ts)) if ts[i] is not None]
        test_scores += [ts[i] for i in xrange(len(ts)) if ts[i] is not None]

    scorer_fprs = []
    scorer_tprs = []
    scorer_precisions = []
    scorer_recalls = []
    for allscores_col in allscores_cols:
        fprs, tprs, _ = roc_curve(test_scores, allscores_col, pos_label=1)
        scorer_fprs.append(fprs)
        scorer_tprs.append(tprs)
        precisions, recalls, _ = precision_recall_curve(test_scores, allscores_col, pos_label=1)
        scorer_precisions.append(precisions)
        scorer_recalls.append(recalls)

    print_auc("PR", batchscore_ids, scorer_recalls, scorer_precisions)
    print_auc("ROC", batchscore_ids, scorer_fprs, scorer_tprs)

    plot_pr(dataset_name, scorer_precisions, scorer_recalls, batchscore_ids)
    plot_roc(dataset_name, scorer_fprs, scorer_tprs, batchscore_ids, .5)
    plt.show(block=False)

    print """\n* To plot another PR curve:
plot_pr(dataset_name, scorer_precisions, scorer_recalls, batchscore_ids, legend='upper right')
plt.show()"""
    print """\n* To plot another ROC curve:
plot_roc(dataset_name, scorer_fprs, scorer_tprs, batchscore_ids, x_max=.5, legend='lower right')
plt.show()"""
    print ""
    import ipdb
    ipdb.set_trace()


def plot_pr(name, scorer_precisions, scorer_recalls, batchscore_ids, legend='upper right'):
    fig = plt.figure()
    y_max = 0
    for ps, rs, id in zip(scorer_precisions, scorer_recalls, batchscore_ids):
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


def plot_roc(name, scorer_fprs, scorer_tprs, batchscore_ids, x_max=1, legend='lower right'):
    fig = plt.figure()

    y_max = 0
    for fprs, tprs, batchscore_id in zip(scorer_fprs, scorer_tprs, batchscore_ids):
        if x_max != 1:
            ind = bisect.bisect(fprs, x_max)
            fprs = fprs[:ind+1]
            tprs = tprs[:ind+1]
            y_max = max(y_max, tprs[-1])
        plt.plot(fprs, tprs, label=batchscore_id)
    plt.plot([0,1],[0,1],'k--')

    if x_max != 1:
        x_rg = [0,x_max]
        y_rg = [0,y_max]
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


AUC_LEVELS = [.1, .5, 1]
def print_auc(name, batchscore_ids, scorer_xs, scorer_ys):
    """
    This doesn't work for PR curves
    """
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
    for batchscore_id, scorer_auc in zip(batchscore_ids, scorers_aucs):
        line = []
        for auc in scorer_auc:
            line.append("%.4f"%auc)
        line.append(batchscore_id)
        print "\t".join(line)
