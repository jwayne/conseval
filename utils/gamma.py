from __future__ import division
import numpy as np
from scipy.special import gammainc, gammaincinv


class DiscreteGammaDistribution(object):

    MAX_RATE = 20

    def __init__(self, alpha, beta, K):
        if alpha <= 0:
            raise Exception("alpha = %f <= 0" % alpha)
        if beta <= 0:
            raise Exception("beta = %f <= 0" % beta)
        if K < 1:
            raise Exception("Too few bins")

        self.alpha = alpha
        self.beta = beta
        self.K = K

        # bin only between [0, MAX_RATE]
        # XXX: this doesn't really seem necessary, since in going from
        # 1 to cum_prob we're really not changing the bin_mids by that much...
        cum_prob = beta * gammainc(alpha, beta*self.MAX_RATE)
        bin_prob = cum_prob / K

        # find expected value of each bin
        # rate4site does this in a weird way that gives slightly different
        # results.  For their implementation, see:
        #   gammaUtilities::the_avarage_r_in_category_between_a_and_b
        #   generalGammaDistribution::fill_mean
        targets = np.arange(bin_prob/2, cum_prob, bin_prob)
        bin_mids = gammaincinv(alpha, targets/beta) / beta

        # Rate of middle of each bin
        self.bin_rates = bin_mids
        # Probability mass of each bin
        self.bin_probs = np.zeros(K) + bin_prob


    def get_rates(self):
        """
        Get the expected values of each bin.
        Note that each bin has the same probability mass, as dictated by
        the discrete gamma model (see Susko 2002)
        """
        return self.bin_rates
