import pylab as pl
from sklearn.metrics import roc_curve, auc


# it = iterator on (alignment, score_tups)

N_BINS = 201

def roc(it, scorer_names):
    """
    build roc curve for each alignment, for each scorer
    """

    # fprs are spaced evenly.
    master_fpr = np.linspace(0,1,N_BINS)
    # List of true positive rates for each scorer.
    master_tprs = []
    for scorer_name in scorer_names:
        master_tprs.append(np.zeros(len(master_fpr)))

    count = 0
    for alignment, score_tups in it:
        scores_cols = zip(*score_tups)
        for master_tpr, scores in zip(tprs, scores_cols):
            # Find this roc curve
            fpr, tpr, thresholds = roc_curve(alignment.testset, scores)
            # Average this roc curve with the master one.
            ind = 0
            for i in xrange(N_BINS):
                while master_fpr[i] < fpr[ind]:
                    ind += 1
                master_tpr[i] += tpr[ind]
        count += 1

    for master_tpr in master_tprs:
        master_tpr /= count

    pl.clf()
    for master_tpr, scorer_name in zip(master_tprs, scorer_names):
        pl.plot(master_fpr, master_tpr, label=scorer_name)
    pl.plot([0,1],[0,1],'k--')
    pl.xlim([0,1])
    pl.ylim([0,1])
    pl.xlabel('False Positive Rate')
    pl.ylabel('True Positive Rate')
    pl.title('ROC curves')
    pl.legend(loc="lower right")
    pl.show()
    #pl.savefig(fname)
