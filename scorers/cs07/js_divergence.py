"""
Jensen-Shannon Divergence (Capra and Singh 07)
Code copyright Tony Capra 2007.
"""
import math
from scorers.cs07.base import Cs07Scorer
from substitution import paramdef_bg_distribution
from utils.bio import weighted_freq_count_pseudocount, PSEUDOCOUNT


class JsDivergence(Cs07Scorer):

    params = Cs07Scorer.params.extend(paramdef_bg_distribution)

    def _score_col(self, col, seq_weights):
        """
        Return the Jensen-Shannon Divergence for the column with the background
        distribution bg_distr.
        """
        distr = self.bg_distribution

        fc = weighted_freq_count_pseudocount(col, seq_weights, PSEUDOCOUNT)

        # if background distrubtion lacks a gap count, remove fc gap count
        if len(distr) == 20:
            new_fc = fc[:-1]
            s = sum(new_fc)
            for i in xrange(len(new_fc)):
                new_fc[i] = new_fc[i] / s
            fc = new_fc

        assert len(fc) == len(distr)

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

        return d
