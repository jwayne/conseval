"""
Code largely by Tony Capra 2007.
"""
from __future__ import division

from params import ParamDef
from scorers.base import Scorer
from utils.bio import get_column, weighted_gap_penalty


class Cs07Scorer(Scorer):

    params = Scorer.params.extend(
        ParamDef('gap_cutoff', .3, float, lambda x: 0<=x<=1,
            help="maximum allowed fraction of gaps per column; columns > this won't be scored"),
        ParamDef('use_gap_penalty', True, bool,
            help="penalize sites with gaps, by multiplying their score by the fraction of non-gap positions in the column."),
        ParamDef('use_seq_weights', True, bool,
            help="weight sequences in alignments to adjust for the overall level of similarity in the alignment."),
    )

    def _score(self, alignment):
        if self.use_seq_weights:
            seq_weights = alignment.get_seq_weights()
        else:
            seq_weights = [1.] * len(alignment.msa)

        scores = []
        for i in xrange(len(alignment.msa[0])):
            col = get_column(i, alignment.msa)
            n_gaps = col.count('-')
            assert n_gaps < len(col)
            if self.gap_cutoff != 1 and n_gaps/len(col) > self.gap_cutoff:
                score = None #want rate for no conservation
            else:
                score = self._score_col(col, seq_weights)
                if self.use_gap_penalty:
                    # vn_entropy has this commented out for some reason
                    score *= weighted_gap_penalty(col, seq_weights)
            scores.append(score)
        return scores


    def _score_col(self, col, seq_weights):
        raise NotImplementedError()
