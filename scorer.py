import math
import numpy as np
import re
import os
import sys

from params import ParamDef, Params
from utils.bio import gap_percentage, get_column, amino_acids, aa_to_index
from substitution import SubstitutionModel, read_sim_matrix, read_bg_distribution


################################################################################
# Get a Scorer by name
################################################################################

def get_scorer(name, **params):
    try:
        scorer_cls = get_scorer_cls(name)
    except (ImportError, AttributeError), e:
        sys.stderr.write("%s: %s is not a valid scoring method.\n" % (e, name))
        return None
    scorer = scorer_cls(**params)
    return scorer

def get_scorer_cls(name):
    scorer_module = __import__('scorers.'+name, fromlist=['x'])
    scorer_clsname = "".join(s.capitalize() for s in name.split('_'))
    return getattr(scorer_module, scorer_clsname)


################################################################################
# Helper classes
################################################################################

class Precache(object):
    pass


################################################################################
# Scorer class
################################################################################

class Scorer(object):

    # Tunable parameters for scoring.  Children inherit all these parameters
    # along with these defaults.  Defaults can be overridden and parameters
    # can be extended, see scorers/mayrose04.py for an example.
    PARAMS = Params(
        #dat matrix file of rate matrix AND bg distribution
        ParamDef('sub_model_file', 'matrix/jtt-dcmut.dat.txt',
            lambda x: os.path.abspath(x)),
        #similarity matrix file, *.bla or *.qij
        ParamDef('sim_matrix_file', 'matrix/blosum62.bla',
            lambda x: os.path.abspath(x)),
        #background distribution file, e.g., swissprot.distribution
        ParamDef('bg_distribution_file', 'matrix/blosum62.distribution',
            lambda x: os.path.abspath(x)),
        #Number of residues on either side included in the window
        ParamDef('window_size', 3, int, lambda x: x>=0),
        #lambda for window heuristic linear combination
        ParamDef('window_lambda', .5, float, lambda x: 0<=x<=1),
        #Do not score columns that contain more than gap cutoff fraction gaps
        ParamDef('gap_cutoff', .3, float, lambda x: 0<=x<=1),
        #Print the z-score (over the alignment) of each column raw score
        #penalize gaps by this amount
        ParamDef('normalize_scores', False, bool),
        #penalize gaps by this amount. The gap penalty used is
        #the score times the fraction of non-gap positions in the column.
        ParamDef('gap_penalty', 1, float),
    )

    # If there are errors.
    DEFAULT_SCORE = 0.

    USE_DAT_MATRIX_AND_DISTRIBUTION = False
    USE_SIM_MATRIX = False
    USE_BG_DISTRIBUTION = False


    def __init__(self, **params):
        s1 = re.sub('(.)([A-Z][a-z]+)', r'\1_\2', type(self).__name__)
        self.name = re.sub('([a-z0-9])([A-Z])', r'\1_\2', s1).lower()

        self.PARAMS.set_params(self, params)

        if self.USE_DAT_MATRIX_AND_DISTRIBUTION and \
                (self.USE_SIM_MATRIX or self.USE_BG_DISTRIBUTION):
            raise Exception()
        # Data file from Kosiol & Goldman 04
        # http://www.ebi.ac.uk/goldman/dayhoff/
        if self.USE_DAT_MATRIX_AND_DISTRIBUTION:
            self.sub_model = SubstitutionModel(self.sub_model_file)
        # Code from Capra & Singh 07
        if self.USE_SIM_MATRIX:
            self.sim_matrix = read_sim_matrix(self.sim_matrix_file)
        # Code from Capra & Singh 07
        if self.USE_BG_DISTRIBUTION:
            self.bg_distribution = read_bg_distribution(self.bg_distribution_file)


    def get_params(self):
        return self.PARAMS.get_params(self)


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
        # Precache computations shared among calls to score_col().
        precache = Precache()
        self._precache(alignment, precache)

        # Main computation.
        scores = []
        for i in range(len(alignment.msa[0])):
            col = get_column(i, alignment.msa)

            if len(col) == len(alignment.msa):
                if self.gap_cutoff == 1 or \
                        gap_percentage(col) <= self.gap_cutoff:
                    scores.append(self.score_col(col, precache))
                else:
                    scores.append(self.DEFAULT_SCORE)
            else:
                sys.stderr.write("Missing sequences in column %d\n" % i)
        #print "Mean score: %f" % np.mean(scores)
        if self.window_size > 0:
            scores = window_score(scores, self.window_size,
                    self.window_lambda)
        if self.normalize_scores:
            scores = calc_z_scores(scores, -999)
        return scores


    def _precache(self, alignment, precache):
        """ 
        Override this function to pre-cache computations shared among calls
        to score_col().  Store any pre-cached computations as attributes of
        `precache`.
        """
        return


    def score_col(self, col, precache):
        """
        @param col:
            the column to be scored.
        @param precache:
            object whose attributes are precache computations specific to this
            alignment.
        """
        raise NotImplementedError()


################################################################################
# Score adjustments
################################################################################

def window_score(scores, window_len, lam=.5):
    """ This function takes a list of scores and a length and transforms them
    so that each position is a weighted average of the surrounding positions.
    Positions with scores less than zero are not changed and are ignored in the
    calculation. Here window_len is interpreted to mean window_len residues on
    either side of the current residue. """

    w_scores = scores[:]

    for i in range(window_len, len(scores) - window_len):
        if scores[i] < 0:
            continue

        sum = 0.
        num_terms = 0.
        for j in range(i - window_len, i + window_len + 1):
            if i != j and scores[j] >= 0:
                num_terms += 1
                sum += scores[j]

        if num_terms > 0:
            w_scores[i] = (1 - lam) * (sum / num_terms) + lam * scores[i]

    return w_scores


def calc_z_scores(scores, score_cutoff):
    """Calculates the z-scores for a set of scores. Scores below
    score_cutoff are not included."""

    average = 0.
    std_dev = 0.
    z_scores = []
    num_scores = 0

    for s in scores:
        if s > score_cutoff:
            average += s
            num_scores += 1
    if num_scores != 0:
        average /= num_scores

    for s in scores:
        if s > score_cutoff:
            std_dev += ((s - average)**2) / num_scores
    std_dev = math.sqrt(std_dev)

    for s in scores:
        if s > score_cutoff and std_dev != 0:
            z_scores.append((s-average)/std_dev)
        else:
            z_scores.append(-1000.0)

    return z_scores
