"""
Relative Entropy (Samudrala and Wang 06)
Code copyright Tony Capra 2007.
"""
import math
from scorers.cs07.base import Cs07Scorer
from substitution import paramdef_bg_distribution
from utils.bio import weighted_freq_count_pseudocount, PSEUDOCOUNT


class RelativeEntropy(Cs07Scorer):

    params = Cs07Scorer.params.extend(paramdef_bg_distribution)

    def _score_col(self, col, seq_weights):
        """
        Calculate the relative entropy of the column distribution with a
        background distribution specified in bg_distr. This is similar to the
        approach proposed in Wang and Samudrala 06.
        """
        distr = self.bg_distribution

        fc = weighted_freq_count_pseudocount(col, seq_weights, PSEUDOCOUNT)

        # remove gap count
        if len(distr) == 20:
            new_fc = fc[:-1]
            s = sum(new_fc)
            for i in xrange(len(new_fc)):
                new_fc[i] = new_fc[i] / s
            fc = new_fc

        assert len(fc) == len(distr)

        d = 0.
        for i in xrange(len(fc)):
            if distr[i] != 0.0:
                d += fc[i] * math.log(fc[i]/distr[i])

        d /= math.log(len(fc))

        return d
