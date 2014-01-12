"""
Empirical Bayes with discretized Gamma prior with the specific
parameters of the gamma determined by a second latent variable, c,
of functionality or non-functionality.
Code by Josh Chen 2013
"""
from __future__ import division
from collections import defaultdict
import numpy as np

from conseval.params import ParamDef
from conseval.scorer import Scorer
from conseval.substitution import paramdef_sub_model
from conseval.utils.bio import get_column
from conseval.utils.gamma import DiscreteGammaDistribution
from scorers.rate4site_eb import compute_subtree_likelihood, precompute_tree_probs


# This is to avoid spurious discrete gamma distributions.
MAX_ALPHA = 40

#XXX: How should I set this?
EM_CONVERGENCE = 100


class Rate4siteFunc(Scorer):

    params = Scorer.params.extend(
#        ParamDef('theta', .05, float, lambda x: 0<x<1,
#            help="proportion of sites that is expected to be functional"),
#        ParamDef('alpha_f', , float, lambda x: 0<x<MAX_ALPHA,
#            help="alpha parameter into the gamma prior on rate r for functional sites."),
#        ParamDef('alpha_nf', 0, float, lambda x: 0<x<MAX_ALPHA,
#            help="alpha parameter into the gamma prior on rate r for non-functional sites."),
        ParamDef('K', 16, int, lambda x: x>1,
            help="number of bins to use for the discrete approximation to the gamma distribution"),
        paramdef_sub_model,
    )


    def __init__(self, **params):
        super(Rate4siteFunc, self).__init__(**params)

        # Precompute this for the initial alpha because it will be used at the
        # start every time an alignment is scored
        theta_1 = .05
        mu_0 = 1.03
        sig2_0 = .5
        mu_1 = .3
        sig2_1 = .06
        self.prior_distr = CombinedDGD(theta_1, mu_0, sig2_0, mu_1, sig2_1, self.K)


    def _score(self, alignment):
        # Silly check; return mean if there aren't enough sequences in the
        # alignment to return reasonable scores.
        if len(alignment.msa) <= 2:
            return [1] * len(alignment.msa[0])

        names_map = dict((name, i) for i, name in enumerate(alignment.names))
        rates, rates_for_est_0, rates_for_est_1, cs, log_marginal = \
                self._estimate_r(alignment, names_map, self.prior_distr)

        # EM for empirical bayes estimate of alpha
        prev_log_marginal = None
        count = 1
        # Perform empirical Bayes estimation of alpha until convergence
        # of the likelihood.
        while not prev_log_marginal or \
                abs(log_marginal-prev_log_marginal) > EM_CONVERGENCE:
            prev_log_marginal = log_marginal
            # M-step
            theta_1 = len(rates_for_est_1) / (len(rates_for_est_1) + len(rates_for_est_0))
            mu_0 = np.mean(rates_for_est_0)
            sig2_0 = np.var(rates_for_est_0)
            mu_1 = np.mean(rates_for_est_1)
            sig2_1 = np.var(rates_for_est_1)
            prior_distr = CombinedDGD(theta_1, mu_0, sig2_0, mu_1, sig2_1, self.K)
            # Compute next E step
            rates, rates_for_est_0, rates_for_est_1, cs, log_marginal = \
                    self._estimate_r(alignment, names_map, prior_distr)
            count += 1

        # negate the rates so higher scores are conserved, just like all the
        # other scorers
        import ipdb
        ipdb.set_trace()
        return [-r for r in rates]


    def _estimate_r(self, alignment, names_map, prior_distr):
        tree = alignment.get_phylotree()

        # Pre-compute the probabilities for every branch and rate.
        # This can be done because the discrete gamma distribution tells us
        # which rates P(rt) will be computed for when scoring columns.
        P_cached = precompute_tree_probs(tree, prior_distr.get_rates(), self.sub_model)

        # E-step
        rates = []
        rates_for_est_0 = []
        rates_for_est_1 = []
        cs = []
        cs_for_est = []
        log_marginal = 0
        for i in xrange(len(alignment.msa[0])):
            col = get_column(i, alignment.msa)
            n_gaps = col.count('-')
            assert n_gaps < len(col)
            if n_gaps == len(col) - 1:
                # Return mean rate.
                rate = 1
                # Return no function.
                c = 0
            else:
                rate, c, marginal = self._estimate_r_col(col, names_map, prior_distr, P_cached, tree)
                if c:
                    rates_for_est_1.append(rate)
                else:
                    rates_for_est_0.append(rate)
                log_marginal += np.log(marginal)
            rates.append(rate)
            cs.append(c)

        return rates, rates_for_est_0, rates_for_est_1, cs, log_marginal


    def _estimate_r_col(self, col, names_map, prior_distr, P_cached, tree):
        """
        Compute this site's rate of evolution r as the expectation of the
        posterior: E[r|X] = \sum_r( P[X|r] P[r] r ) / \sum_r( P[X|r] P[r] ).

        Assume a fixed alpha until I get it working (leave the inference of the
        hyperparameter until later).
        """
        root = tree.root

        # Numerator \sum_r,c( r*P(X,r,c) )
        top_r = 0
        # Numerator \sum_r,c( c*P(X,r,c) )
        top_c = 0
        # Denominator \sum_r,c( P(X,r,c) )
        marginal = 0
        for p, c, r in prior_distr.get_probs_rates():
            # P(X|r)
            likelihood = compute_subtree_likelihood(root, r, col, names_map, P_cached)
            # likelihood None only if column is all gaps
            assert likelihood is not None
            # likelihood might be returned as a 1x1 matrix
            likelihood = float(likelihood)
            # P(X,r,c).
            joint = likelihood * p
            top_r += joint * r
            if c:
                top_c += joint
            marginal += joint
        exp_r = top_r / marginal
        exp_c = int(top_c / marginal * 2)
        return exp_r, exp_c, marginal


class CombinedDGD(object):

    def __init__(self, theta, alpha_0, beta_0, alpha_1, beta_1, K):
        if not (0 <= theta <= 1):
            raise Exception("theta = %f not in (0,1)" % theta)
        self.theta = theta
        self.dgd_0 = DiscreteGammaDistribution(alpha_0, beta_0, K)
        self.dgd_1 = DiscreteGammaDistribution(alpha_1, beta_1, K)
        self.K = K
        self.probs_rates = [(1-self.theta, 0, r) for r in self.dgd_0.get_rates()] + \
                           [(self.theta, 1, r) for r in self.dgd_1.get_rates()]

    def get_rates(self):
        """
        Get rates without their probabilities.  Different rates will have different
        probabilities.  This should be used only for precomputing P_cached.
        """
        return [r for _, _, r in self.probs_rates]

    def get_probs_rates(self):
        """
        Get rates with their probabilities, and their c.
        """
        return self.probs_rates
