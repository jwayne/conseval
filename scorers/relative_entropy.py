"""
Relative Entropy (Samudrala and Wang 06)
Code copyright Tony Capra 2007.
"""
import math
from scorer import Scorer
from utils.bio import weighted_freq_count_pseudocount, weighted_gap_penalty, PSEUDOCOUNT


class RelativeEntropy(Scorer):

    USE_BG_DISTRIBUTION = True

    def _precache(self, alignment, precache):
        precache.seq_weights = alignment.get_seq_weights()

    def score_col(self, col, precache):
        """
        Calculate the relative entropy of the column distribution with a
        background distribution specified in bg_distr. This is similar to the
        approach proposed in Wang and Samudrala 06.
        """
        distr = self.bg_distribution[:]

        fc = weighted_freq_count_pseudocount(col, precache.seq_weights, PSEUDOCOUNT)

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

        if self.gap_penalty == 1:
            return d * weighted_gap_penalty(col, precache.seq_weights)
        else:
            return d


