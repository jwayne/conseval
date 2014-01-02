"""
Property Relative Entropy (Williamson 95)
Code copyright Tony Capra 2007.
"""
import math
from scorers.cs07.base import Cs07Scorer
from substitution import paramdef_bg_distribution
from utils.bio import aa_to_index, weighted_freq_count_pseudocount, PSEUDOCOUNT


class PropertyRelativeEntropy(Cs07Scorer):

    params = Cs07Scorer.params.extend(paramdef_bg_distribution)

    property_partition = [['V','L', 'I','M'], ['F','W','Y'], ['S','T'], ['N','Q'], ['H','K','R'], ['D','E'], ['A','G'], ['P'], ['C']]

    def _score_col(self, col, seq_weights):
        """
        Calculate the relative entropy of a column col relative to a
        partition of the amino acids. Similar to Williamson '95.  See shannon_entropy()
        for more general info.
        """
        prop_bg_freq = []
        bg_distr = self.bg_distribution
        if len(bg_distr) == len(self.property_partition):
            prop_bg_freq = bg_distr
        else:
            prop_bg_freq = [0.248, 0.092, 0.114, 0.075, 0.132, 0.111, 0.161, 0.043, 0.024, 0.000] # from BL62

        fc = weighted_freq_count_pseudocount(col, seq_weights, PSEUDOCOUNT)

        # sum the aa frequencies to get the property frequencies
        prop_fc = [0.] * len(self.property_partition)
        for p in xrange(len(self.property_partition)):
            for aa in self.property_partition[p]:
                prop_fc[p] += fc[aa_to_index[aa]]

        d = 0.
        for i in xrange(len(prop_fc)):
            if prop_fc[i] != 0 and prop_bg_freq[i] != 0:
                d += prop_fc[i] * math.log(prop_fc[i] / prop_bg_freq[i], 2)

        return d
