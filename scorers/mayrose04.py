"""
Empirical Bayes with Gamma prior, binned into 4 categories (Mayrose et al 04)
See http://www.tau.ac.il/~itaymay/cp/rate4site.html for the original version.
Code by Josh Chen 2013
"""
from __future__ import division
from collections import defaultdict
import numpy as np

from scorer import Scorer, ParamDef
from utils.bio import aa_to_index
from utils.gamma import DiscreteGammaDistribution


class Mayrose04(Scorer):

    PARAMS = Scorer.PARAMS.override({
        'window_size': 0,
        'gap_cutoff': 1,
    }).extend(
        ParamDef('alpha', 1, float, lambda x: x>0,
            help="alpha parameters into the gamma prior"),
        ParamDef('n_gamma_bins', 16, int, lambda x: x>0,
            help="number of bins to use for the discrete approximation to the gamma distribution"),
    )

    USE_DAT_MATRIX_AND_DISTRIBUTION = True

    def __init__(self, **params):
        super(Mayrose04, self).__init__(**params)
        self.beta = self.alpha
        self.dgd = DiscreteGammaDistribution(
                self.alpha, self.beta, self.n_gamma_bins)
        self.DEFAULT_SCORE = self.alpha / self.beta


    def _precache(self, alignment, precache):
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

        precache.names_map = dict((name, i) for i, name in enumerate(alignment.names))
        precache.P_cached = P_cached
        precache.tree = tree


    def score_col(self, col, precache):
        """
        Compute this site's rate of evolution r as the expectation of the
        posterior: E[r|X] = \sum_r( P[X|r] P[r] r ) / \sum_r( P[X|r] P[r] ).

        Assume a fixed alpha until I get it working (leave the inference of the
        hyperparameter until later).
        """
        names_map = precache.names_map
        P_cached = precache.P_cached
        tree = precache.tree
        root = tree.clade

        # Compute numerator and denominator separately
        top = 0
        bot = 0
        for rate, prior in self.dgd.get_categories():
            likelihood = self._compute_subtree_likelihood(root, rate, col, names_map, P_cached)
            # likelihood None only if column is all gaps
            assert likelihood is not None
            joint = float(likelihood) * prior
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
            if aa == '-':
                # Ignored gapped columns
                return None
            P_node_given_parent = P_cached[rate][node][:,aa_to_index[aa]]
            return P_node_given_parent
        else:
            child_P_cols = []
            for child in node.clades:
                child_P_col = self._compute_subtree_likelihood(child, rate, col, names_map, P_cached)
                if child_P_col is not None:
                    child_P_cols.append(child_P_col)
            if not child_P_cols:
                return None
            P_children_given_node = np.matrix(np.prod(child_P_cols, axis=0))
            P_subtree_given_parent = P_cached[rate][node] * P_children_given_node
            return P_subtree_given_parent
