from __future__ import division
import bisect
import numpy as np
import pylab as pl
from sklearn.metrics import roc_curve

from utils.stats import zscore


AUC_LEVELS = [.1, .5, 1]


# it = iterator on (alignment, score_tups)

def roc(it, scorer_names, out_dir):
    """
    build roc curve for each alignment, for each scorer
    """
    scorers_scores = []
    for i in scorer_names:
        scorers_scores.append([])
    test_scores = []

    # Just aggregate all scores across all data files.  It isn't much memory anyway.
    count = 0
    for alignment, score_tups in it:
        scores_cols = zip(*score_tups)
        for scorer_scores, scores in zip(scorers_scores, scores_cols):
            scorer_scores += scores
        test_scores += alignment.testset

    scorer_fprs = []
    scorer_tprs = []
    for scorer_scores in scorers_scores:
        fprs, tprs, _ = roc_curve(test_scores, scorer_scores, pos_label=1)
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

    pl.clf()
    for fprs, tprs, scorer_name in zip(scorer_fprs, scorer_tprs, scorer_names):
        pl.plot(fprs, tprs, label=scorer_name)
    pl.plot([0,1],[0,1],'k--')
    pl.xlim([0,1])
    pl.ylim([0,1])
    pl.xlabel('False Positive Rate')
    pl.ylabel('True Positive Rate')
    pl.title('ROC curves')
    pl.legend(loc="lower right")
    pl.show()
    #pl.savefig(fname)
