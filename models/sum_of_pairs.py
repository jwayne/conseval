from models.scorer import Scorer
from utils import *


################################################################################
# Mutation Weighted Pairwise Match
################################################################################

class SumOfPairs(Scorer):

    USE_SIM_MATRIX = True

    def score_col(self, col, seq_weights, gap_penalty=1):
        """ Sum the similarity matrix values for all pairs in the column.
        This method is similar to those proposed in Valdar 02."""

        sum = 0.
        max_sum = 0.

        for i in range(len(col)):
            for j in range(i):
                if col[i] != '-' and col[j] != '-':
                    max_sum += seq_weights[i] * seq_weights[j]
                    sum += seq_weights[i] * seq_weights[j] * \
                            self.sim_matrix[aa_to_index[col[i]]][aa_to_index[col[j]]]

        if max_sum != 0:
            sum /= max_sum
        else:
            sum = 0.

        if gap_penalty == 1:
            return sum * weighted_gap_penalty(col, seq_weights)
        else:
            return sum
