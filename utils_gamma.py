from __future__ import division
from scipy.special import gammainc


class DiscreteGammaDistribution(object):

    def __init__(self, alpha, beta, n_cats):
        if alpha <= 0:
            raise Exception("alpha = %f <= 0" % alpha)
        if beta <= 0:
            raise Exception("beta = %f <= 0" % beta)
        if n_cats < 1:
            raise Exception("Too few categories")
        n_iters = 0
        max_err = 1e-5

        lo = 0
        cat_ubounds = [0] * n_cats
        for cat in xrange(n_cats):
            target = (cat+1) / n_cats
            hi = lo + .1

            # find suitable lo and hi points, via exponential intervals
            while True:
                n_iters += 1
                curr = gammainc(alpha, hi)
                if curr < target:
                    hi += (hi - lo)
                else:
                    lo += (hi - lo) / 2
                    break
            # find suitable mid point, via binary search
            while True:
                n_iters += 1
                mid = (lo + hi) / 2
                curr = gammainc(alpha, mid)
                err = abs(curr - target)
                if err < max_err:
                    break
                else:
                    if curr > target:
                        hi = mid
                    else:
                        lo = mid
            # add as upper boundary after adjusting by beta
            cat_ubounds[cat] = mid / beta
            # set lo for next round
            lo = mid

        cat_rates = [0] * n_cats
        cat_probs = [0] * n_cats
        lo = 0
        for cat in xrange(n_cats):
            # rate4site represents the bin's rate as the average rate in the bin
            cat_rates[cat] = (cat_ubounds[cat]*beta + lo*beta) / 2
            # rate4site represents the bin's prior probability as the probability mass
            # of the bin
            #XXX: rate4site uses alpha+1??
            #XXX: rate4site multiplies P[r_i] by alpha/beta*k
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

