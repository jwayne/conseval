import math
from models.scorer import Scorer
from utils import *


################################################################################
# Jensen-Shannon Divergence
################################################################################

class JsDivergence(Scorer):

    USE_BG_DISTRIBUTION = True

    def score_col(self, col, seq_weights, gap_penalty=1):
        """Return the Jensen-Shannon Divergence for the column with the background
        distribution bg_distr."""

        distr = self.bg_distribution[:]

        fc = weighted_freq_count_pseudocount(col, seq_weights, PSEUDOCOUNT)

        # if background distrubtion lacks a gap count, remove fc gap count
        if len(distr) == 20:
            new_fc = fc[:-1]
            s = sum(new_fc)
            for i in range(len(new_fc)):
                new_fc[i] = new_fc[i] / s
            fc = new_fc

        if len(fc) != len(distr): return -1

        # make r distriubtion
        r = []
        for i in range(len(fc)):
            r.append(.5 * fc[i] + .5 * distr[i])

        d = 0.
        for i in range(len(fc)):
            if r[i] != 0.0:
                if fc[i] == 0.0:
                    d += distr[i] * math.log(distr[i]/r[i], 2)
                elif distr[i] == 0.0:
                    d += fc[i] * math.log(fc[i]/r[i], 2)
                else:
                    d += fc[i] * math.log(fc[i]/r[i], 2) + distr[i] * math.log(distr[i]/r[i], 2)

        # d /= 2 * math.log(len(fc))
        d /= 2

        if gap_penalty == 1:
            return d * weighted_gap_penalty(col, seq_weights)
        else:
            return d
