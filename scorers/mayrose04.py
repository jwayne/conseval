"""
Empirical Bayes with Gamma prior, binned into 4 categories (Mayrose et al 04)
See http://www.tau.ac.il/~itaymay/cp/rate4site.html for the original version.
Code by Josh Chen 2013
"""
from __future__ import division
from collections import defaultdict
import numpy as np
from scorer import Scorer
from utils import aa_to_index
from utils_gamma import DiscreteGammaDistribution


class Mayrose04(Scorer):

    PARAM_OVERRIDES = {
        'window_size': 0,
        'gap_cutoff': 1,
    }

    USE_DAT_MATRIX_AND_DISTRIBUTION = True

    ALPHA = 1
    BETA = 1
    N_GAMMA_CATEGORIES = 16

    def __init__(self, **params):
        super(Mayrose04, self).__init__(**params)
        self.dgd = DiscreteGammaDistribution(
                self.ALPHA, self.BETA, self.N_GAMMA_CATEGORIES)
        self.DEFAULT_SCORE = self.ALPHA / self.BETA


    def _precache(self, alignment, precached):
        P_cached = defaultdict(dict)
        tree = alignment.get_phylotree()

        terminals = set()
        root = tree.clade
        bfs = [root]
        for node in bfs:
            for rate, _ in self.dgd.get_categories():
                if node == root:
                    P_cached[rate][node] = self.sub_model.freqs
                else:
                    t = node.branch_length
                    P = self.sub_model.calc_P(rate*t)
                    P_cached[rate][node] = P
            if not node.is_terminal():
                bfs += node.clades
            else:
                # Sanity check.
                if not node.name:
                    raise Exception("Node has no name")
                if node.name in terminals:
                    raise Exception("Duplicate node name: %s" % node.name)
                terminals.add(node.name)

        precached.names_map = dict((name, i) for i, name in enumerate(alignment.names))
        precached.P_cached = P_cached
        precached.tree = tree


    def score_col(self, col, precached):
        """
        Compute this site's rate of evolution r as the expectation of the
        posterior: E[r|X] = \sum_r( P[X|r] P[r] r ) / \sum_r( P[X|r] P[r] ).

        Assume a fixed alpha until I get it working (leave the inference of the
        hyperparameter until later).
        """
        if '-' in col:
            # Return average rate if there exists a gap
            # TODO
            return self.DEFAULT_SCORE

        names_map = precached.names_map
        P_cached = precached.P_cached
        tree = precached.tree
        root = tree.clade

        # Compute numerator and denominator separately
        top = 0
        bot = 0
        for rate, prior in self.dgd.get_categories():
            likelihood = self._compute_subtree_likelihood(root, rate, col, names_map, P_cached)
            joint = likelihood * prior
            top += joint * rate
            bot += joint
        expectation = top / bot
        return expectation


    def _compute_subtree_likelihood(self, node, rate, col, names_map, P_cached):
        """
        Helper function to compute likelihood via postorder traversal.
        """
        if node.is_terminal():
            aa = col[names_map[node.name]]
            P_node_given_parent = P_cached[rate][node][:,aa_to_index[aa]]
            return P_node_given_parent
        else:
            child_P_cols = [self._compute_subtree_likelihood(child, rate, col, names_map, P_cached)
                    for child in node.clades]
            P_children_given_node = np.matrix(np.prod(child_P_cols, axis=0))
            P_subtree_given_parent = P_cached[rate][node] * P_children_given_node
            return P_subtree_given_parent
