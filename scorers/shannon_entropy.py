"""
Shannon Entropy
Code copyright Tony Capra 2007.
"""
import math
from scorer import Scorer
from utils.bio import weighted_freq_count_pseudocount, weighted_gap_penalty, PSEUDOCOUNT


class ShannonEntropy(Scorer):

    def _precache(self, alignment, precache):
        precache.seq_weights = alignment.get_seq_weights()

    def score_col(self, col, precache):
        """
        Calculates the Shannon entropy of the column col.
        The entropy will be between zero and one because of its base. See p.13 of
        Valdar 02 for details. The information score 1 - h is returned for the sake
        of consistency with other scores.
        """
        fc = weighted_freq_count_pseudocount(col, precache.seq_weights, PSEUDOCOUNT)

        h = 0.
        for i in range(len(fc)):
            if fc[i] != 0:
                h += fc[i] * math.log(fc[i])

        #h /= math.log(len(fc))
        h /= math.log(min(len(fc), len(col)))

        inf_score = 1 - (-1 * h)

        if self.gap_penalty == 1:
            return inf_score * weighted_gap_penalty(col, precache.seq_weights)
        else:
            return inf_score
