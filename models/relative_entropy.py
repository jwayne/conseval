import math
from models.scorer import Scorer
from utils import *


################################################################################
# Relative Entropy
################################################################################

class RelativeEntropy(Scorer):

    USE_BG_DISTRIBUTION = True

    def score_col(self, col, seq_weights, gap_penalty=1):
        """Calculate the relative entropy of the column distribution with a
        background distribution specified in bg_distr. This is similar to the
        approach proposed in Wang and Samudrala 06."""

        distr = self.bg_distribution[:]

        fc = weighted_freq_count_pseudocount(col, seq_weights, PSEUDOCOUNT)

        # remove gap count
        if len(distr) == 20:
            new_fc = fc[:-1]
            s = sum(new_fc)
            for i in range(len(new_fc)):
                new_fc[i] = new_fc[i] / s
            fc = new_fc

        if len(fc) != len(distr): return -1

        d = 0.
        for i in range(len(fc)):
            if distr[i] != 0.0:
                d += fc[i] * math.log(fc[i]/distr[i])

        d /= math.log(len(fc))

        if gap_penalty == 1:
            return d * weighted_gap_penalty(col, seq_weights)
        else:
            return d


