from __future__ import division
import numpy as np
from scipy.special import gammainc, gammaincinv


class DiscreteGammaDistribution(object):

    # This is actually needed so that gammaincinv doesn't give inf
    MAX_RATE = 20

    def __init__(self, alpha, beta, K):
        if alpha <= 0:
            raise Exception("alpha = %f <= 0" % alpha)
        if beta <= 0:
            raise Exception("beta = %f <= 0" % beta)
        if K < 1:
            raise Exception("Too few categories")

        # find upper boundaries of each category
        max_prob = gammainc(alpha, self.MAX_RATE)
        bin_prob = max_prob / K
        targets = np.arange(bin_prob, max_prob+bin_prob/2, bin_prob)
        #XXX: not sure why we don't divide targets / beta
        cat_ubounds = gammaincinv(alpha, targets) / beta

        cat_lbounds = np.zeros(K)
        cat_lbounds[1:] = cat_ubounds[:-1]
        tmp = gammainc(alpha+1, cat_ubounds * beta) - gammainc(alpha+1, cat_lbounds * beta)
        cat_rates = tmp * alpha / beta * K

        # Rate of middle of each category
        self.cat_rates = cat_rates
        # Probability mass of each category
        self.cat_probs = np.zeros(K) + bin_prob


    def get_rates(self):
        return self.cat_rates
