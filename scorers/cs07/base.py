"""
Code largely by Tony Capra 2007.
"""
from __future__ import division

from params import ParamDef
from scorers.base import Scorer
from utils.bio import get_column, weighted_gap_penalty


class Cs07Scorer(Scorer):

    params = Scorer.params.with_defaults({
        'window_size': 3,
    }).extend(
        ParamDef('gap_cutoff', .3, float, lambda x: 0<=x<=1,
            help="maximum allowed fraction of gaps per column; columns > this won't be scored"),
        ParamDef('use_gap_penalty', True, bool,
            help="penalize sites with gaps, by multiplying their score by the fraction of non-gap positions in the column."),
        ParamDef('use_seq_weights', True, bool,
            help="weight sequences in alignments to adjust for the overall level of similarity in the alignment."),
    )

    # If column has more gaps than gap cutoff, give this score. This score
    # should be a score representing "no conservation", since if a column
    # has a lot of gaps we assume that it isn't conserved.
    SCORE_OVER_GAP_CUTOFF = None


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
                score = self.SCORE_OVER_GAP_CUTOFF
            else:
                score = self._score_col(col, seq_weights)
                if self.use_gap_penalty:
                    # vn_entropy has this commented out for some reason
                    score *= weighted_gap_penalty(col, seq_weights)
            scores.append(score)
        return scores


    def _score_col(self, col, seq_weights):
        raise NotImplementedError()
