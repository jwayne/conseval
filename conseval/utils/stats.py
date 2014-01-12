import numpy as np


def zscore(x, filter=0):
    x = np.array(x)
    avg = np.mean(x)
    stdev = np.std(x)
    if not stdev:
        return x
    z_scores = (x - avg) / stdev
    if filter:
        for i in xrange(len(z_scores)):
            if abs(z_scores[i]) > filter:
                z_scores[i] = filter
    return z_scores
