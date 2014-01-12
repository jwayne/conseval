from __future__ import division
import numpy as np
from scipy.special import gammainc, gammaincinv


class DiscreteGammaDistribution(object):
    """
    Discrete approximation to the gamma distribution, using equal probability
    mass bins.
    """

    # This is actually needed so that gammaincinv doesn't give inf
    MAX_RATE = 20

    def __init__(self, alpha, beta, K):
        """
        Note alpha/beta = mean, alpha/beta^2 = variance
        @param alpha:
            shape parameter
        @param beta:
            inverse scale parameter
        @param K:
            number of bins
        """
        if alpha <= 0:
            raise Exception("alpha = %f <= 0" % alpha)
        if beta <= 0:
            raise Exception("beta = %f <= 0" % beta)
        if K < 1:
            raise Exception("Num bins = %d < 1" % K)

        # find upper boundaries of each bin
        max_prob = gammainc(alpha, self.MAX_RATE)
        bin_prob = max_prob / K
        targets = np.arange(bin_prob, max_prob+bin_prob/2, bin_prob)
        #XXX: not sure why we don't divide targets / beta
        bin_ubounds = gammaincinv(alpha, targets) / beta

        bin_lbounds = np.zeros(K)
        bin_lbounds[1:] = bin_ubounds[:-1]
        tmp = gammainc(alpha+1, bin_ubounds * beta) - gammainc(alpha+1, bin_lbounds * beta)
        bin_rates = tmp * alpha / beta * K

        # Rate of middle of each bin
        self.bin_rates = bin_rates
        # Probability mass of each bin
        self.bin_probs = np.zeros(K) + bin_prob


    def get_rates(self):
        return self.bin_rates
