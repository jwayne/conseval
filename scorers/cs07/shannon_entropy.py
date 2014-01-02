"""
Shannon Entropy
Code copyright Tony Capra 2007.
"""
import math
from scorers.cs07.base import Cs07Scorer
from utils.bio import weighted_freq_count_pseudocount, PSEUDOCOUNT


class ShannonEntropy(Cs07Scorer):

    def _score_col(self, col, seq_weights):
        """
        Calculates the Shannon entropy of the column col.
        The entropy will be between zero and one because of its base. See p.13 of
        Valdar 02 for details. The information score 1 - h is returned for the sake
        of consistency with other scores.
        """
        fc = weighted_freq_count_pseudocount(col, seq_weights, PSEUDOCOUNT)

        h = 0.
        for i in xrange(len(fc)):
            if fc[i] != 0:
                h += fc[i] * math.log(fc[i])

        #h /= math.log(len(fc))
        h /= math.log(min(len(fc), len(col)))

        inf_score = 1 - (-1 * h)

        return inf_score
