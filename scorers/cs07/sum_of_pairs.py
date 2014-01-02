"""
Mutation Weighted Pairwise Match
Code copyright Tony Capra 2007.
"""
from scorers.cs07.base import Cs07Scorer
from substitution import paramdef_sim_matrix
from utils.bio import aa_to_index


class SumOfPairs(Cs07Scorer):

    params = Cs07Scorer.params.extend(paramdef_sim_matrix)

    def _score_col(self, col, seq_weights):
        """
        Sum the similarity matrix values for all pairs in the column.
        This method is similar to those proposed in Valdar 02.
        """
        curr_sum = 0.
        max_sum = 0.
        for i in xrange(len(col)):
            for j in xrange(i):
                if col[i] != '-' and col[j] != '-':
                    max_sum += seq_weights[i] * seq_weights[j]
                    curr_sum += seq_weights[i] * seq_weights[j] * \
                            self.sim_matrix[aa_to_index[col[i]]][aa_to_index[col[j]]]

        if max_sum != 0:
            curr_sum /= max_sum
        else:
            curr_sum = 0.

        return curr_sum
