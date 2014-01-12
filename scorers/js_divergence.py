"""
Jensen-Shannon Divergence (Capra and Singh 07)
Code by Josh Chen 2013.  Idea from Tony Capra 2007.
"""
import math
import numpy as np
from conseval.params import ParamDef
from conseval.scorer import Scorer
from conseval.substitution import paramdef_bg_distribution_optional
from conseval.utils.bio import weighted_freq_count_pseudocount, PSEUDOCOUNT, amino_acids, get_column, weighted_gap_penalty


class JsDivergence(Scorer):

    params = Scorer.params.with_defaults({
        'window_size': 2,
    }).extend(
        ParamDef('gap_cutoff', .3, float, lambda x: 0<=x<=1,
            help="maximum allowed fraction of gaps per column; columns > this won't be scored"),
        ParamDef('use_gap_penalty', True, bool,
            help="penalize sites with gaps, by multiplying their score by the fraction of non-gap positions in the column."),
        ParamDef('use_seq_weights', True, bool,
            help="weight sequences in alignments to adjust for the overall level of similarity in the alignment."),
        ParamDef('lambda_prior', .1, float, lambda x: 0<=x<=1,
            help="prior weight lambda_prior in the Jensen-Shannon divergence"),
        paramdef_bg_distribution_optional
    )

    SCORE_OVER_GAP_CUTOFF = 0


    def _score(self, alignment):
        if self.use_seq_weights:
            seq_weights = alignment.get_seq_weights()
        else:
            seq_weights = [1.] * len(alignment.msa)

        if self.bg_distribution is None:
            # Estimate bg distribution from this alignment
            q = weighted_freq_count_pseudocount((aa for seq in alignment.msa for aa in seq),
                    seq_weights, PSEUDOCOUNT)
        else:
            q = self.bg_distribution

        scores = []
        for i in xrange(len(alignment.msa[0])):
            col = get_column(i, alignment.msa)
            n_gaps = col.count('-')
            assert n_gaps < len(col)
            if self.gap_cutoff != 1 and n_gaps/len(col) > self.gap_cutoff:
                score = self.SCORE_OVER_GAP_CUTOFF
            else:
                score = self._score_col(col, seq_weights, q)
                if self.use_gap_penalty:
                    # vn_entropy has this commented out for some reason
                    score *= weighted_gap_penalty(col, seq_weights)
            scores.append(score)
        return scores


    def _score_col(self, col, seq_weights, q):
        """
        Return the Jensen-Shannon Divergence for the column with the background
        distribution q.
        """
        lamb1 = self.lambda_prior
        lamb2 = 1-self.lambda_prior

        # get frequency distribution
        with_gap = (len(q) == 21)
        pc = weighted_freq_count_pseudocount(col, seq_weights, PSEUDOCOUNT, with_gap)
        assert len(pc) == len(q)

        # make r distriubtion
        r = lamb1*pc + lamb2*q

        # sum relative entropies
        d1 = lamb1 * sum(pc[i] * math.log(pc[i]/r[i], 2) for i in xrange(len(pc)) if pc[i])
        d2 = lamb2 * sum(q[i] * math.log(q[i]/r[i], 2) for i in xrange(len(pc)) if pc[i])

        return d1+d2
