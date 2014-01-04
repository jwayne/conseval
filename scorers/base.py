import numpy as np
import os
import sys
import time

from params import ParamDef, Params, WithParams
from utils.stats import zscore


################################################################################
# Get a Scorer by name
################################################################################

def get_scorer(name, **params):
    """
    Get the appropriate Scorer object in the module scorers.`name`, initialized
    with 'params'.
    """
    scorer_cls = get_scorer_cls(name)
    scorer = scorer_cls(**params)
    return scorer

def get_scorer_cls(name):
    """
    Get the appropriate Scorer class in the module scorers.`name`.
      >>> cls = get_scorer_cls('caprasingh07.js_divergence')
      >>> cls.__module__
      scorers.caprasingh07.js_divergence
      >>> cls.__name__
      JsDivergence
    """
    try:
        scorer_module = __import__('scorers.'+name, fromlist=['x'])
        scorer_clsname = "".join(s.capitalize() for s in name.split('.')[-1].split('_'))
        scorer_cls = getattr(scorer_module, scorer_clsname)
    except (ImportError, AttributeError), e:
        raise ImportError("%s: %s is not a valid scorer." % (e, name))
    return scorer_cls


################################################################################
# Scorer class
################################################################################

class Scorer(WithParams):
    """
    Base class for any scorer.  Define subclasses that override the _score
    method in scorers/[scorer_name.py] with class name ScorerName.
    """

    # Tunable parameters for scoring.  Children inherit all these parameters
    # along with these defaults.  Defaults can be overridden and parameters
    # can be extended, see scorers/mayrose04.py for an example.
    params = Params(
        ParamDef('window_size', 0, int, lambda x: x>=0,
            help="Number of residues on either side included in the window"),
        ParamDef('window_lambda', .5, float, lambda x: 0<=x<=1,
            help="lambda for window heuristic linear combination. Meaningful only if window_size != 0."),
        ParamDef('normalize', False, bool,
            help="return z-scores (over the alignment) of each column, instead of original scores"),
    )


    def __init__(self, **params):
        super(Scorer, self).__init__(**params)
        self.name = ".".join(type(self).__module__.split('.')[1:])


    def score(self, alignment):
        """
        Score each site in the first sequence of `alignment`.  Performs computations
        that are not specific to any site, and calls score_col() to perform the
        site-secific computations.

        Additional global computations can be performed by overriding _precache(),
        see below.

        @param alignment:
            Alignment object
        @return:
            List of scores for each site
        """
        t0 = time.time()

        # Main computation.
        scores = self._score(alignment)

        if self.window_size:
            scores = window_score(scores, self.window_size,
                    self.window_lambda)
        if self.normalize:
            scores = list(zscore(scores))

        dt = time.time() - t0 #len(alignment.msa), len(alignment.msa[0])
        return scores


    def _score(self, alignment):
        """
        Called by _score(..).  Override to define scoring method for entire
        alignment.

        @param alignment:
            Alignment object
        @return:
            List of scores for each site
        """
        raise NotImplementedError()




################################################################################
# Score adjustments
################################################################################

def window_score(scores, window_len, lam=.5):
    """
    This function takes a list of scores and a length and transforms them
    so that each position is a weighted average of the surrounding positions.
    Positions with scores less than zero are not changed and are ignored in the
    calculation. Here window_len is interpreted to mean window_len residues on
    either side of the current residue.
    
    Code by Tony Capra 2007.
    """
    w_scores = scores[:]

    for i in xrange(window_len, len(scores) - window_len):
        if scores[i] < 0:
            continue

        curr_sum = 0.
        num_terms = 0.
        for j in xrange(i - window_len, i + window_len + 1):
            if i != j and scores[j] >= 0:
                num_terms += 1
                curr_sum += scores[j]

        if num_terms > 0:
            w_scores[i] = (1 - lam) * (curr_sum / num_terms) + lam * scores[i]

    return w_scores
