"""
Mutation Weighted Pairwise Match
Code copyright Tony Capra 2007.
"""
from scorers.cs07.base import Cs07Scorer
from conseval.substitution import paramdef_sim_matrix
from conseval.utils.bio import aa_to_index


class SumOfPairs(Cs07Scorer):

    # Normalize because the score doesn't go from 0 to 1...
    params = Cs07Scorer.params.with_defaults({
        'normalize': True
    }).extend(paramdef_sim_matrix)

    # TODO: determine what score is representative of non-conserved scores,
    # for SCORE_OVER_GAP_CUTOFF


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

        if max_sum:
            curr_sum /= max_sum
        else:
            curr_sum = 0.

        return curr_sum
