from __future__ import division
from scipy.special import gammainc, gammaincinv


class DiscreteGammaDistribution(object):

    MAX_RATE = 20

    def __init__(self, alpha, beta, n_cats):
        if alpha <= 0:
            raise Exception("alpha = %f <= 0" % alpha)
        if beta <= 0:
            raise Exception("beta = %f <= 0" % beta)
        if n_cats < 1:
            raise Exception("Too few categories")
        n_iters = 0

        # find upper boundaries of each category
        MAX_PROB = gammainc(alpha, self.MAX_RATE)
        cat_ubounds = []
        for cat in xrange(n_cats):
            target = (cat+1) / n_cats * MAX_PROB
            cat_ubounds.append(gammaincinv(alpha, target) / beta)

        cat_rates = [0] * n_cats
        cat_probs = [0] * n_cats
        lo = 0
        for cat in xrange(n_cats):
            # rate4site represents the bin's rate as the average rate in the bin
            cat_rates[cat] = (cat_ubounds[cat]*beta + lo*beta) / 2
            # rate4site represents the bin's prior probability as the probability mass
            # of the bin
            #XXX: why does rate4site use alpha+1?
            #XXX: why does rate4site multiply P[r_i] by alpha/beta*n_cats?
            #XXX: I get uniform cat_probs without the above 2 mods, which makes sense.
            # I don't get why rate4site does something different then.
            cat_probs[cat] = 1/n_cats
            cat_probs[cat] = (gammainc(alpha+1, cat_ubounds[cat]*beta) - \
                    gammainc(alpha+1, lo*beta)) * alpha / beta * n_cats
            # set lo for next round
            lo = cat_ubounds[cat]

        # Rate of upper bound of each category
        self.cat_ubounds = cat_ubounds
        # Rate of middle of each category
        self.cat_rates = cat_rates
        # Probability mass of each category
        self.cat_probs = cat_probs
        # Number of iterations it took to find the upper bounds
        self.n_iters = n_iters


    def get_categories(self):
        for rate, prob in zip(self.cat_rates, self.cat_probs):
            yield rate, prob


    def __str__(self):
        return "ubounds %s\nrates %s\nprobs %s\n%d iterations" % (
            self.cat_ubounds,
            self.cat_rates,
            self.cat_probs,
            self.n_iters)
