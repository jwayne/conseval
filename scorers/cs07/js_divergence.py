"""
Jensen-Shannon Divergence (Capra and Singh 07)
Code by Josh Chen 2013.  Idea from Tony Capra 2007.
"""
import math
from params import ParamDef
from scorers.cs07.base import Cs07Scorer
from substitution import paramdef_bg_distribution
from utils.bio import weighted_freq_count_pseudocount, PSEUDOCOUNT


class JsDivergence(Cs07Scorer):

    params = Cs07Scorer.params.extend(
        ParamDef('prior_weight', .5, float, lambda x: 0<=x<=1,
            help="prior weight lambda in the Jensen-Shannon divergence"),
        paramdef_bg_distribution
    )

    SCORE_OVER_GAP_CUTOFF = 0


    def _score_col(self, col, seq_weights):
        """
        Return the Jensen-Shannon Divergence for the column with the background
        distribution q.
        """
        q = self.bg_distribution
        lamb1 = self.prior_weight
        lamb2 = 1-self.prior_weight

        # get frequency distribution
        with_gap = (len(q) == 21)
        pc = weighted_freq_count_pseudocount(col, seq_weights, PSEUDOCOUNT, with_gap)
        assert len(pc) == len(q)

        # make r distriubtion
        r = lamb1*pc + lamb2*q

        # sum relative entropies
        d1 = lamb1 * sum(pc[i] * math.log(pc[i]/r[i], 2) for i in xrange(len(pc)) if pc[i])
        d2 = lamb2 * sum(q[i] * math.log(q[i]/r[i], 2) for i in xrange(len(pc)) if pc[i])

        return d1+d2
