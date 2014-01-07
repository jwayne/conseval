from __future__ import division
import bisect
import numpy as np
import matplotlib.pyplot as plt
from sklearn.metrics import roc_curve

from evaluate import get_batchscores, get_out_dir
from conseval.utils.stats import zscore


AUC_LEVELS = [.1, .5, 1]


# it = iterator on (alignment, score_tups)

def roc(dataset_name, *scorer_ids):
    """
    build roc curve for each alignment, for each scorer
    """
    scores_cols = [[] for i in scorer_ids]
    test_scores = []

    # Just aggregate all scores across all data files.  It isn't much memory anyway.
    for alignment, scores_tup in get_batchscores(dataset_name, scorer_ids):
        ts = alignment.testset
        for scores_col, scores in zip(scores_cols, scores_tup):
            scores_col += [scores[i] for i in xrange(len(ts)) if ts[i] is not None]
        test_scores += [s for s in ts if s is not None]

    scorer_fprs = []
    scorer_tprs = []
    for scores_col in scores_cols:
        fprs, tprs, _ = roc_curve(test_scores, scores_col, pos_label=1)
        scorer_fprs.append(fprs)
        scorer_tprs.append(tprs)

    # Compute AUCs for each scorer
    scorers_aucs = []
    for fprs, tprs in zip(scorer_fprs, scorer_tprs):
        scorer_aucs = []
        for auc_level in AUC_LEVELS:
            ind = bisect.bisect(fprs, auc_level)
            auc = np.trapz(tprs[:ind], fprs[:ind])
            scorer_aucs.append(auc)
        scorers_aucs.append(scorer_aucs)

    print ""
    print "Scorer\t%s" % "\t".join("AUC_{%d}"%int(aucl*100) for aucl in AUC_LEVELS)
    for scorer_id, scorer_auc in zip(scorer_ids, scorers_aucs):
        line = [scorer_id]
        for auc in scorer_auc:
            line.append("%.4f"%auc)
        print "\t".join(line)

    plt.figure()
    for fprs, tprs, scorer_id in zip(scorer_fprs, scorer_tprs, scorer_ids):
        plt.plot(fprs, tprs, label=scorer_id)
    plt.plot([0,1],[0,1],'k--')
    plt.xlim([0,1])
    plt.ylim([0,1])
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.title('ROC curves')
    plt.legend(loc="lower right")

    plt.figure()
    fpr_max = 0.05
    tpr_max = 0
    for fprs, tprs, scorer_id in zip(scorer_fprs, scorer_tprs, scorer_ids):
        ind = bisect.bisect(fprs, fpr_max)
        plt.plot(fprs[:ind+1], tprs[:ind+1], label=scorer_id)
        tpr_max = max(tpr_max, np.max(tprs[:ind+1]))
    plt.plot([0,1],[0,1],'k--')
    plt.xlim([0,fpr_max])
    plt.ylim([0,tpr_max])
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.title('ROC curves')
    plt.legend(loc="lower right")

    plt.show()
