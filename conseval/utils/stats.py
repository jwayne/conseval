import numpy as np


def zscore(x):
    x = np.array(x)
    avg = np.mean(x)
    stdev = np.std(x)
    if not stdev:
        return x
    z_scores = (x - avg) / stdev
    return z_scores
