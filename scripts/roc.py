from __future__ import division
import bisect
import numpy as np
import pylab as pl
from sklearn.metrics import roc_curve

from evaluate import get_batchscores, get_out_dir
from conseval.utils.stats import zscore


AUC_LEVELS = [.1, .5, 1]


# it = iterator on (alignment, score_tups)

def roc(dataset_name, scorer_ids, id=None):
    """
    build roc curve for each alignment, for each scorer
    """
    out_dir = get_out_dir(id)

    scores_cols = [[] for i in scorer_ids]
    test_scores = []

    # Just aggregate all scores across all data files.  It isn't much memory anyway.
    for alignment, scores_tup in get_batchscores(dataset_name, scorer_ids):
        for scores_col, scores in zip(scores_cols, scores_tup):
            scores_col += scores
        test_scores += alignment.testset

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

    pl.clf()
    for fprs, tprs, scorer_name in zip(scorer_fprs, scorer_tprs, scorer_ids):
        pl.plot(fprs, tprs, label=scorer_name)
    pl.plot([0,1],[0,1],'k--')
    pl.xlim([0,1])
    pl.ylim([0,1])
    pl.xlabel('False Positive Rate')
    pl.ylabel('True Positive Rate')
    pl.title('ROC curves')
    pl.legend(loc="lower right")
    pl.show()
    #l.savefig(fname)
