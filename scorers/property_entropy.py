"""
Property Entropy (Mirny and Shakhnovich 95, Valdar and Thornton 01)
Code copyright Tony Capra 2007.
"""
import math
from scorer import Scorer
from utils import *


class PropertyEntropy(Scorer):

    def score_col(self, col, seq_weights, gap_penalty=1, alignment=None):
        """Calculate the entropy of a column col relative to a partition of the
        amino acids. Similar to Mirny '99."""

        # Mirny and Shakn. '99
        property_partition = [['A','V','L','I','M','C'], ['F','W','Y','H'], ['S','T','N','Q'], ['K','R'], ['D', 'E'], ['G', 'P'], ['-']]

        # Williamson '95
        # property_partition = [['V','L', 'I','M'], ['F','W','Y'], ['S','T'], ['N','Q'], ['H','K','R'], ['D','E'], ['A','G'], ['P'], ['C'], ['-']]

        fc = weighted_freq_count_pseudocount(col, seq_weights, PSEUDOCOUNT)

        # sum the aa frequencies to get the property frequencies
        prop_fc = [0.] * len(property_partition)
        for p in range(len(property_partition)):
            for aa in property_partition[p]:
                prop_fc[p] += fc[aa_to_index[aa]]

        h = 0.
        for i in range(len(prop_fc)):
            if prop_fc[i] != 0:
                h += prop_fc[i] * math.log(prop_fc[i])

        h /= math.log(min(len(property_partition), len(col)))

        inf_score = 1 - (-1 * h)

        if gap_penalty == 1:
            return inf_score * weighted_gap_penalty(col, seq_weights)
        else:
            return inf_score
