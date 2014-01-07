"""
Empirical Bayes with discretized Gamma prior (Mayrose et al 04)
See http://www.tau.ac.il/~itaymay/cp/rate4site.html for the original version.
Code by Josh Chen 2013
"""
from __future__ import division
from collections import defaultdict
import numpy as np

from conseval.params import ParamDef
from conseval.scorer import Scorer
from conseval.substitution import paramdef_sub_model
from conseval.utils.bio import aa_to_index, get_column
from conseval.utils.gamma import DiscreteGammaDistribution


class Rate4siteEb(Scorer):

    params = Scorer.params.extend(
        ParamDef('alpha', 0, float, lambda x: x>=0,
            help="alpha parameter into the gamma prior on rate r. Var(r) = 1/alpha. If alpha=0 (default), then the empirical Bayesian estimate of alpha is used."),
        ParamDef('K', 16, int, lambda x: x>1,
            help="number of bins to use for the discrete approximation to the gamma distribution"),
        paramdef_sub_model,
    )

    #XXX: How should I set this?
    EM_CONVERGENCE = 100


    def __init__(self, **params):
        super(Rate4siteEb, self).__init__(**params)

        if self.alpha:
            # Fixed alpha
            alpha = self.alpha
        else:
            # Initial estimate of alpha for empirical Bayes estimation
            alpha = 1
        # Precompute this for the initial alpha because it will be used at the
        # start every time an alignment is scored
        self.dgd = DiscreteGammaDistribution(alpha, alpha, self.K)


    def _score(self, alignment):
        names_map = dict((name, i) for i, name in enumerate(alignment.names))
        rates, alpha_est, log_marginal = self._score_fixed_alpha(
                alignment, names_map, self.dgd)

        # Fixed alpha
        if self.alpha:
            return rates
        # Bad alpha_est (due to all sites having too many gaps)
        if alpha_est is None:
            return rates

        # EM for empirical bayes estimate of alpha
        prev_log_marginal = None
        while not prev_log_marginal or \
                    abs(log_marginal-prev_log_marginal) > self.EM_CONVERGENCE:
            prev_log_marginal = log_marginal
            dgd = DiscreteGammaDistribution(alpha_est, alpha_est, self.K)
            rates, alpha_est, log_marginal = self._score_fixed_alpha(
                    alignment, names_map, dgd)
        return rates


    def _score_fixed_alpha(self, alignment, names_map, dgd):
        tree = alignment.get_phylotree()

        # Pre-compute the probabilities for every branch and rate.
        # This can be done because the discrete gamma distribution tells us
        # which rates P(rt) will be computed for when scoring columns.
        P_cached = defaultdict(dict)
        root = tree.clade
        bfs = [root]
        for node in bfs:
            for rate in dgd.get_rates():
                if node is root:
                    P_cached[rate][node] = self.sub_model.freqs
                else:
                    t = node.branch_length
                    P = self.sub_model.calc_P(rate*t)
                    P_cached[rate][node] = P

            if not node.is_terminal():
                bfs += node.clades

        # E-step
        rates = []
        rates_for_est = []
        log_marginal = 0
        for i in xrange(len(alignment.msa[0])):
            col = get_column(i, alignment.msa)
            n_gaps = col.count('-')
            assert n_gaps < len(col)
            if n_gaps == len(col) - 1:
                # Return mean rate.
                rate = 1
            else:
                rate, marginal = self._score_col(col, names_map, dgd, P_cached, tree)
                rates_for_est.append(rate)
                log_marginal += np.log(marginal)
            rates.append(rate)

        # M-step
        if len(rates_for_est) <= 1:
            alpha_est = None
        else:
            alpha_est = 1 / np.var(rates_for_est)

        return rates, alpha_est, log_marginal


    def _score_col(self, col, names_map, dgd, P_cached, tree):
        """
        Compute this site's rate of evolution r as the expectation of the
        posterior: E[r|X] = \sum_r( P[X|r] P[r] r ) / \sum_r( P[X|r] P[r] ).

        Assume a fixed alpha until I get it working (leave the inference of the
        hyperparameter until later).
        """
        root = tree.clade

        # Numerator \sum_r( r*P(X,r) )
        top = 0
        # Denominator \sum_r( P(X,r) )
        bot = 0
        for rate in dgd.get_rates():
            # P(X|r)
            likelihood = self._compute_subtree_likelihood(root, rate, col, names_map, P_cached)
            # likelihood None only if column is all gaps
            assert likelihood is not None
            # likelihood might be returned as a 1x1 matrix
            likelihood = float(likelihood)
            # P(X,r).  Since in the discrete gamma model, the probability of each
            # bin is the same, we don't multiply by the prior.
            joint = likelihood
            top += joint * rate
            bot += joint
        expectation = top / bot
        return expectation, bot


    def _compute_subtree_likelihood(self, node, rate, col, names_map, P_cached):
        """
        Helper function to compute likelihood via postorder traversal.
        """
        if node.is_terminal():
            aa = col[names_map[node.name]]
            if aa is '-':
                # Ignored nodes set to gaps
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
