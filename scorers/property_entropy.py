"""
Property Entropy (Mirny and Shakhnovich 95, Valdar and Thornton 01)
Code copyright Tony Capra 2007.
"""
import math
from scorer import Scorer
from utils.bio import (aa_to_index,
                   weighted_freq_count_pseudocount, weighted_gap_penalty, PSEUDOCOUNT)


class PropertyEntropy(Scorer):

    # Mirny and Shakn. '99
    property_partition = [['A','V','L','I','M','C'], ['F','W','Y','H'], ['S','T','N','Q'], ['K','R'], ['D', 'E'], ['G', 'P'], ['-']]

    def _precache(self, alignment, precache):
        precache.seq_weights = alignment.get_seq_weights()

    def score_col(self, col, precache):
        """
        Calculate the entropy of a column col relative to a partition of the
        amino acids. Similar to Mirny '99.
        """
        fc = weighted_freq_count_pseudocount(col, precache.seq_weights, PSEUDOCOUNT)

        # sum the aa frequencies to get the property frequencies
        prop_fc = [0.] * len(self.property_partition)
        for p in range(len(self.property_partition)):
            for aa in self.property_partition[p]:
                prop_fc[p] += fc[aa_to_index[aa]]

        h = 0.
        for i in range(len(prop_fc)):
            if prop_fc[i] != 0:
                h += prop_fc[i] * math.log(prop_fc[i])

        h /= math.log(min(len(self.property_partition), len(col)))

        inf_score = 1 - (-1 * h)

        if self.gap_penalty == 1:
            return inf_score * weighted_gap_penalty(col, precache.seq_weights)
        else:
            return inf_score
