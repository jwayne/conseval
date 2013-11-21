import math
from scorer import Scorer
from utils import *


################################################################################
# Property Relative Entropy
################################################################################

class PropertyRelativeEntropy(Scorer):

    USE_BG_DISTRIBUTION = True

    def score_col(self, col, seq_weights, gap_penalty=1):
        """Calculate the relative entropy of a column col relative to a
        partition of the amino acids. Similar to Williamson '95.  See shannon_entropy()
        for more general info. """

        # Mirny and Shakn. '99
        #property_partition = [['A','V','L','I','M','C'], ['F','W','Y','H'], ['S','T','N','Q'], ['K','R'], ['D', 'E'], ['G', 'P'], ['-']]

        # Williamson '95
        property_partition = [['V','L', 'I','M'], ['F','W','Y'], ['S','T'], ['N','Q'], ['H','K','R'], ['D','E'], ['A','G'], ['P'], ['C']]

        prop_bg_freq = []
        bg_distr = self.bg_distribution
        if len(bg_distr) == len(property_partition):
            prop_bg_freq = bg_distr
        else:
            prop_bg_freq = [0.248, 0.092, 0.114, 0.075, 0.132, 0.111, 0.161, 0.043, 0.024, 0.000] # from BL62

        #fc = weighted_freq_count_ignore_gaps(col, seq_weights)
        fc = weighted_freq_count_pseudocount(col, seq_weights, PSEUDOCOUNT)

        # sum the aa frequencies to get the property frequencies
        prop_fc = [0.] * len(property_partition)
        for p in range(len(property_partition)):
            for aa in property_partition[p]:
                prop_fc[p] += fc[aa_to_index[aa]]

        d = 0.
        for i in range(len(prop_fc)):
            if prop_fc[i] != 0 and prop_bg_freq[i] != 0:
                d += prop_fc[i] * math.log(prop_fc[i] / prop_bg_freq[i], 2)


        if gap_penalty == 1:
            return d * weighted_gap_penalty(col, seq_weights)
        else:
            return d
