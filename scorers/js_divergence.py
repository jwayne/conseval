"""
Jensen-Shannon Divergence (Capra and Singh 07)
Code copyright Tony Capra 2007.
"""
import math
from scorer import Scorer
from utils.bio import weighted_freq_count_pseudocount, weighted_gap_penalty, PSEUDOCOUNT


class JsDivergence(Scorer):

    USE_BG_DISTRIBUTION = True

    def _precache(self, alignment, precache):
        precache.seq_weights = alignment.get_seq_weights()

    def score_col(self, col, precache):
        """
        Return the Jensen-Shannon Divergence for the column with the background
        distribution bg_distr.
        """
        distr = self.bg_distribution[:]

        fc = weighted_freq_count_pseudocount(col, precache.seq_weights, PSEUDOCOUNT)

        # if background distrubtion lacks a gap count, remove fc gap count
        if len(distr) == 20:
            new_fc = fc[:-1]
            s = sum(new_fc)
            for i in xrange(len(new_fc)):
                new_fc[i] = new_fc[i] / s
            fc = new_fc

        if len(fc) != len(distr): return -1

        # make r distriubtion
        r = []
        for i in xrange(len(fc)):
            r.append(.5 * fc[i] + .5 * distr[i])

        d = 0.
        for i in xrange(len(fc)):
            if r[i] != 0.0:
                if fc[i] == 0.0:
                    d += distr[i] * math.log(distr[i]/r[i], 2)
                elif distr[i] == 0.0:
                    d += fc[i] * math.log(fc[i]/r[i], 2)
                else:
                    d += fc[i] * math.log(fc[i]/r[i], 2) + distr[i] * math.log(distr[i]/r[i], 2)

        # d /= 2 * math.log(len(fc))
        d /= 2

        if self.gap_penalty == 1:
            return d * weighted_gap_penalty(col, precache.seq_weights)
        else:
            return d
