from __future__ import division
import matplotlib.pyplot as plt
import numpy as np
import random
from evaluate import get_batchscores


def compare_scores(dataset_name, *batchscore_ids):
    """
    Compare scores from 2 scoring methods in a scatter plot.
    """
    if len(batchscore_ids) != 2:
        raise Exception("Need 2 batchscore runs to compare")
    id1, id2 = batchscore_ids
    pos_pairs = []
    neg_pairs = []

    for alignment, scores_cols in get_batchscores(dataset_name, batchscore_ids):
        ts = alignment.testset
        scores_pairs = zip(*scores_cols)
        pos_pairs += (scores_pairs[i] for i in xrange(len(scores_pairs)) if ts[i])
        neg_pairs += (scores_pairs[i] for i in xrange(len(scores_pairs)) if not ts[i])

    plt.figure()
    x1,y1 = zip(*random.sample(neg_pairs,min(len(neg_pairs),10000)))
    plt.scatter(x1,y1,color='red')
    x2,y2 = zip(*random.sample(pos_pairs,min(len(pos_pairs),1000)))
    plt.scatter(x2,y2,color='green')
    plt.xlabel(id1)
    plt.ylabel(id2)
    plt.xlim([min(min(x1),min(x2)), max(max(x1),max(x2))])
    plt.ylim([min(min(y1),min(y2)), max(max(y1),max(y2))])
    plt.show()
