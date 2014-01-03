"""
Property Entropy (Mirny and Shakhnovich 95, Valdar and Thornton 01)
Code copyright Tony Capra 2007.
"""
import math
from scorers.cs07.base import Cs07Scorer
from utils.bio import aa_to_index, weighted_freq_count_pseudocount, PSEUDOCOUNT


class PropertyShannonEntropy(Cs07Scorer):

    # Mirny and Shakn. '99
    property_partition = [['A','V','L','I','M','C'], ['F','W','Y','H'], ['S','T','N','Q'], ['K','R'], ['D', 'E'], ['G', 'P'], ['-']]

    SCORE_OVER_GAP_CUTOFF = 0


    def _score_col(self, col, seq_weights):
        """
        Calculate the entropy of a column col relative to a partition of the
        amino acids. Similar to Mirny '99.
        """
        fc = weighted_freq_count_pseudocount(col, seq_weights, PSEUDOCOUNT)

        # sum the aa frequencies to get the property frequencies
        prop_fc = [0.] * len(self.property_partition)
        for p in range(len(self.property_partition)):
            for aa in self.property_partition[p]:
                prop_fc[p] += fc[aa_to_index[aa]]

        h = 0.
        for pfc_i in prop_fc:
            if pfc_i:
                h -= pfc_i * math.log(pfc_i)

        # Convert score so that it's between 0 and 1.
        # Recall that shannon entropy is between 0 and log(number of values with nonzero freq)
        # XXX: Why involve len(col) if we have a pseudocount?
        h /= math.log(min(len(self.property_partition), len(col)))

        # Convert score so that 1 is conserved, and 0 is not.
        return 1 - h
