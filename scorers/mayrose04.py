from __future__ import division
import numpy as np
from scorer import Scorer
from utils import aa_to_index
from utils_gamma import DiscreteGammaDistribution


################################################################################
# Empirical Bayes with Gamma prior, binned into 4 categories
################################################################################

class Mayrose04(Scorer):

    USE_DAT_MATRIX_AND_DISTRIBUTION = True
    ALPHA = 1
    BETA = 1
    N_GAMMA_CATEGORIES = 16

    def __init__(self, *args, **kwargs):
        super(Mayrose04, self).__init__(*args, **kwargs)
        self.dgd = DiscreteGammaDistribution(
                self.ALPHA, self.BETA, self.N_GAMMA_CATEGORIES)


    def score(self, alignment, **kwargs):
        alignment.names_map = dict((name, i) for i, name in enumerate(alignment.names))
        tree = alignment.get_phylotree()
        tree_probs = []
#            S = self.sim_matrix
#            PI = self.bg_distribution
#            Q = S*(PI as matrix) ?
#            Q = S ?
#            A = PI^(1/2)*S*PI(1/2)
#            X, Lamb = np.eigen(A)
        all_terminals = []
        for rate, prior in self.dgd.get_categories():
            root = tree.clade
            likelihood_arr = self.bg_distribution
            bfs = [(child, child.branch_length, likelihood_arr) for child in root.clades]
            terminals = {}
            for node, t, likelihood_arr in bfs:
#                P = e^(Q*t) = PI^(-1/2) * X*e^(Lamb*t)*X.T * PI^(1/2)
                likelihood_arr = P * likelihood_arr
                if node.is_terminal():
                    if not node.name:
                        raise Exception("No node name")
                    if node.name in terminals:
                        raise Exception("Duplicate node name: %s" % node.name)
                    terminals[node.name] = likelihood_arr
            all_terminals.append((rate, prior, terminals))
        # TODO: don't save this as instance variable, but pass into score_col instead.
        self.all_terminals = all_terminals
        return super(Mayrose04, self).score(alignment, **kwargs)


    # TODO: make alignment mandatory
    def score_col(self, col, seq_weights, gap_penalty=1, alignment=None):
        """
        Compute this site's rate of evolution r as the expectation of the
        posterior: E[r|X] = \sum_r( P[X|r] P[r] r ) / \sum_r( P[X|r] P[r] ).

        Assume a fixed alpha until I get it working (leave the inference of the
        hyperparameter until later).
        """
        names_map = alignment.names_map

        # Compute numerator and denominator separately
        top = 0
        bot = 0
        for rate, prior, likelihood_arr in self.all_terminals:
            likelihood = 1
            for aa, seq_name in zip(col, names_map):
                aa_index = aa_to_index[aa]
                likelihood *= likelihood_arr[seq_name][aa_index]
            joint = likelihood * prior
            top += joint * rate
            bot += joint
        expectation = top / bot
        return expectation
