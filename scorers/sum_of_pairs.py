"""
Mutation Weighted Pairwise Match
Code copyright Tony Capra 2007.
"""
from scorer import Scorer
from utils.bio import aa_to_index, weighted_gap_penalty


class SumOfPairs(Scorer):

    USE_SIM_MATRIX = True

    def _precache(self, alignment, precache):
        precache.seq_weights = alignment.get_seq_weights()

    def score_col(self, col, precache):
        """
        Sum the similarity matrix values for all pairs in the column.
        This method is similar to those proposed in Valdar 02.
        """
        seq_weights = self.seq_weights

        sum = 0.
        max_sum = 0.
        for i in xrange(len(col)):
            for j in xrange(i):
                if col[i] != '-' and col[j] != '-':
                    max_sum += seq_weights[i] * seq_weights[j]
                    sum += seq_weights[i] * seq_weights[j] * \
                            self.sim_matrix[aa_to_index[col[i]]][aa_to_index[col[j]]]

        if max_sum != 0:
            sum /= max_sum
        else:
            sum = 0.

        if self.gap_penalty == 1:
            return sum * weighted_gap_penalty(col, seq_weights)
        else:
            return sum
