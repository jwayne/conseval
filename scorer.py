import math
import numpy as np
import os
import sys
from utils import gap_percentage, get_column, amino_acids, aa_to_index
from substitution import SubstitutionModel, read_sim_matrix, read_bg_distribution


################################################################################
# Get a Scorer by name
################################################################################

DEFAULT_SCORER = 'js_divergence'
def parse_scorer_names(scorer_names):
    if not scorer_names:
        return [DEFAULT_SCORER]
    return scorer_names

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

class ParamDef(object):

    def __init__(self, name, default, parse_fxn=None, check_fxn=None, help=None):
        self.name = name
        self.parse_fxn = parse_fxn
        self.check_fxn = check_fxn
        self.default = self.parse(default)

    def parse(self, *args):
        if args:
            val = args[0]
            if self.parse_fxn:
                val = self.parse_fxn(val)
            if self.check_fxn and not self.check_fxn(val):
                raise ValueError("Input %r is unacceptable for parameter %s"
                        % (args[0], self.name))
            return val
        return self.default


class Params(object):

    def __init__(self, param_defs, *overrides):
        """
        Set attributes of this object to be the parameters in `param_defs`,
        overriding default values with the values in `overrides`, preferring
        earlier-specified overrides to later-specified ones.

        Check that all parameters specified in the `overrides` dicts are
        present in `param_defs`.
        """
        param_def_keys = set([pd.name for pd in param_defs])
        for i, override in enumerate(overrides):
            for k in override:
                if k not in param_def_keys:
                    raise AttributeError("Input param %s not allowed for this Scorer" % k)
        for param_def in param_defs:
            k = param_def.name
            found = False
            for override in overrides:
                if k in override:
                    setattr(self, k, param_def.parse(override[k]))
                    found = True
                    break
            if not found:
                setattr(self, k, param_def.parse())

    def __str__(self):
        return str(self.__dict__)

    def __repr__(self):
        return "Params(%s)" % self


class Precache(object):
    pass


################################################################################
# Scorer class
################################################################################

class Scorer(object):

    # Tunable parameters for scoring.
    PARAMS = (
        #dat matrix file of rate matrix AND bg distribution
        ParamDef('dat_file', 'matrix/jtt-dcmut.dat.txt'),
        #similarity matrix file, *.bla or *.qij
        ParamDef('sim_matrix_file', 'matrix/blosum62.bla'),
        #background distribution file, e.g., swissprot.distribution
        ParamDef('bg_distribution_file', 'matrix/blosum62.distribution'),
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
    PARAM_OVERRIDES = {}

    # If there are errors.
    DEFAULT_SCORE = 0

    USE_DAT_MATRIX_AND_DISTRIBUTION = False
    USE_SIM_MATRIX = False
    USE_BG_DISTRIBUTION = False


    def __init__(self, **params):
        self.params = Params(self.PARAMS, params, self.PARAM_OVERRIDES)

        if self.USE_DAT_MATRIX_AND_DISTRIBUTION and \
                (self.USE_SIM_MATRIX or self.USE_BG_DISTRIBUTION):
            raise Exception()
        # Data file from Kosiol & Goldman 04
        # http://www.ebi.ac.uk/goldman/dayhoff/
        if self.USE_DAT_MATRIX_AND_DISTRIBUTION:
            self.sub_model = SubstitutionModel(self.params.dat_file)
        # Code from Capra & Singh 07
        if self.USE_SIM_MATRIX:
            self.sim_matrix = read_sim_matrix(self.params.sim_matrix_file)
        # Code from Capra & Singh 07
        if self.USE_BG_DISTRIBUTION:
            self.bg_distribution = read_bg_distribution(self.params.bg_distribution_file)


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
        precached = Precache()
        self._precache(alignment, precached)

        # Main computation.
        scores = []
        for i in range(len(alignment.msa[0])):
            col = get_column(i, alignment.msa)

            if len(col) == len(alignment.msa):
                if self.params.gap_cutoff == 1 or \
                        gap_percentage(col) <= self.params.gap_cutoff:
                    scores.append(self.score_col(col, precached))
                else:
                    scores.append(self.DEFAULT_SCORE)
            else:
                sys.stderr.write("Missing sequences in column %d\n" % i)
        #print "Mean score: %f" % np.mean(scores)
        if self.params.window_size > 0:
            scores = window_score(scores, self.params.window_size,
                    self.params.window_lambda)
        if self.params.normalize_scores:
            scores = calc_z_scores(scores, -999)
        return scores


    def _precache(self, alignment, precached):
        """ 
        Override this function to pre-cache computations shared among calls
        to score_col().  Store any pre-cached computations as attributes of
        `precached`.
        """
        return


    def score_col(self, col, precached):
        """
        @param col:
            the column to be scored.
        @param precached:
            object whose attributes are precached computations specific to this
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
