import numpy as np
from scorer import Scorer
from utils_gamma import DiscreteGammaDistribution


################################################################################
# Empirical Bayes with Gamma prior, binned into 4 categories
################################################################################

class Mayrose04(Scorer):

    USE_SIM_MATRIX = True
    USE_BG_DISTRIBUTION = True
    SKIP_ADJUSTMENTS = True
    ALPHA = 1
    BETA = 1
    N_GAMMA_CATEGORIES = 16

    def __init__(self, *args, **kwargs):
        super(Mayrose04, self).__init__(*args, **kwargs)
        self.dgd = DiscreteGammaDistribution(
                self.ALPHA, self.BETA, self.N_GAMMA_CATEGORIES)


    def score(self, alignment, **kwargs):
        alignment.names_map = dict((name, i) for i, name in enumerate(alignment.names))
        return super(Mayrose04, self).score(alignment, **kwargs)


    def score_col(self, col, seq_weights, gap_penalty=1, alignment=None):
        """
        Compute this site's rate of evolution r as the expectation of the
        posterior: E[r|X] = \sum_r( P[X|r] P[r] r ) / \sum_r( P[X|r] P[r] ).

        Assume a fixed alpha until I get it working (leave the inference of the
        hyperparameter until later).
        """
        import ipdb
        ipdb.set_trace()
        tree = alignment.get_phylotree()
        names_map = alignment.names_map

        # Compute numerator and denominator separately
        top = 0
        bot = 0
        # Array of joints. For each node, for each possible assignment to that node.
        joints = np.zeros(len(names_map))
        for rate, prior in self.dgd.get_categories():
            # joint = likelihood * prior; we have prior so compute likelihood
            likelihood = compute_likelihood(col, names_map)
            node = tree.root

            alignment.names, alignment.msa

            tree.root
            joint = likelihood * prior
            top += joint * rate
            bot += joint
        expectation = top / bot

