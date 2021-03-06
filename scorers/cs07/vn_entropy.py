"""
Von Neumann Entropy (Caffrey et al 04)
Code copyright Tony Capra 2007.
"""
import math
import numpy as np
from scorers.cs07.base import Cs07Scorer
from conseval.substitution import paramdef_sim_matrix
from conseval.utils.bio import aa_to_index


class VnEntropy(Cs07Scorer):

    params = Cs07Scorer.params.extend(paramdef_sim_matrix)

    SCORE_OVER_GAP_CUTOFF = 0


    def _score_col(self, col, seq_weights):
        """
        Calculate the von Neuman Entropy as described in Caffrey et al. 04.
        This code was adapted from the implementation found in the PFAAT project
        available on SourceForge.
        """

        aa_counts = [0.] * 20
        for aa in col:
            if aa != '-': aa_counts[aa_to_index[aa]] += 1

        dm_aas = []
        for i in xrange(len(aa_counts)):
            if aa_counts[i] != 0:
                dm_aas.append(i)
        dm_size = len(dm_aas)

        # No non-gap amino acids in col
        if dm_size == 0: return 0.0

        row_i = 0
        col_i = 0
        dm = np.zeros((dm_size, dm_size), np.float32)
        for i in xrange(dm_size):
            row_i = dm_aas[i]
            for j in xrange(dm_size):
                col_i = dm_aas[j]
                dm[i][j] = aa_counts[row_i] * self.sim_matrix[row_i][col_i]

        ev = np.linalg.eig(dm)[0].real

        temp = 0.
        for e in ev:
            temp += e

        if temp != 0:
            for i in xrange(len(ev)):
                ev[i] = ev[i] / temp

        vne = 0.0
        for e in ev:
            if e > (10**-10):
                vne -= e * math.log(e) / math.log(20)

        # Presumably this makes it so that 1 is conserved, and 0 is not.
        return 1 - vne
