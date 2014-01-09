"""
Shannon Entropy
Code copyright Tony Capra 2007.
"""
import math
from scorers.cs07.base import Cs07Scorer
from conseval.utils.bio import weighted_freq_count_pseudocount, PSEUDOCOUNT


class ShannonEntropy(Cs07Scorer):

    SCORE_OVER_GAP_CUTOFF = 0


    def _score_col(self, col, seq_weights):
        """
        Calculates the Shannon entropy of the column col.
        The entropy will be between zero and one because of its base. See p.13 of
        Valdar 02 for details. The information score 1 - h is returned for the sake
        of consistency with other scores.
        """
        fc = weighted_freq_count_pseudocount(col, seq_weights, PSEUDOCOUNT)

        h = 0.
        for fc_i in fc:
            if fc_i:
                h -= fc_i * math.log(fc_i)

        # Convert score so that it's between 0 and 1.
        # Recall that shannon entropy is between 0 and log(number of values with nonzero freq)
        # XXX: Why involve len(col) if we have a pseudocount?
        h /= math.log(len(fc))#math.log(min(len(fc), len(col)))

        # Convert score so that 1 is conserved, and 0 is not.
        return 1 - h
