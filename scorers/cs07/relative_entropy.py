"""
Relative Entropy (Samudrala and Wang 06)
Code copyright Tony Capra 2007.
"""
import numpy as np
from scorers.cs07.base import Cs07Scorer
from substitution import paramdef_bg_distribution
from utils.bio import weighted_freq_count_pseudocount, PSEUDOCOUNT


class RelativeEntropy(Cs07Scorer):

    params = Cs07Scorer.params.extend(paramdef_bg_distribution)

    SCORE_OVER_GAP_CUTOFF = 0


    def _score_col(self, col, seq_weights):
        """
        Calculate the relative entropy of the column distribution with a
        background distribution specified in bg_distr. This is similar to the
        approach proposed in Wang and Samudrala 06.
        """
        q = self.bg_distribution

        with_gap = (len(q) == 21)
        fc = weighted_freq_count_pseudocount(col, seq_weights, PSEUDOCOUNT, with_gap)
        assert len(fc) == len(q)

        d = np.sum(fc * np.log(fc/q))

        # Convert score so that it's between 0 and 1.
        # XXX: why is relative entropy assumed to be bounded?
        d /= np.log(len(fc))

        return d
