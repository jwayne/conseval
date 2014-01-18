import numpy as np


################################################################################
# Score adjustments
################################################################################

def norm_scores(x, filter=0):
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


def window_scores(scores, window_size, lam=.5):
    """
    This function takes a list of scores and a length and transforms them
    so that each position is a weighted average of the surrounding positions.
    Positions with scores less than zero are not changed and are ignored in the
    calculation. Here window_size is interpreted to mean window_size residues on
    either side of the current residue.
    
    Code by Tony Capra 2007.
    """
    w_scores = scores[:]
    for i in xrange(window_size, len(scores)-window_size):
        if scores[i] is None:
            continue
        curr_sum = 0.
        num_terms = 0
        for j in xrange(i-window_size, i+window_size+1):
            if i != j and scores[j] is not None:
                curr_sum += scores[j]
                num_terms += 1
        if num_terms:
            w_scores[i] = (1-lam) * (curr_sum/num_terms) + lam * scores[i]
    return w_scores
