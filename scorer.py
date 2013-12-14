import math
import numpy as np
import os
import sys
from utils import gap_percentage, get_column, amino_acids, aa_to_index


# BLOSUM62 background distribution
blosum_background_distr = [0.078, 0.051, 0.041, 0.052, 0.024, 0.034, 0.059, 0.083, 0.025, 0.062, 0.092, 0.056, 0.024, 0.044, 0.043, 0.059, 0.055, 0.014, 0.034, 0.072]

current_dir = os.path.dirname(os.path.realpath(__file__))


################################################################################
# Get a Scorer by name
################################################################################

#defaults
def get_scorer(name, dat_file=None, s_matrix_file=None, bg_distribution_file=None):
    try:
        scorer_cls = get_scorer_cls(name)
    except (ImportError, AttributeError), e:
        sys.stderr.write("%s: %s is not a valid scoring method.\n" % (e, name))
        return None
    scorer = scorer_cls(dat_file, s_matrix_file, bg_distribution_file)
    return scorer

def get_scorer_cls(name):
    scorer_module = __import__('scorers.'+name, fromlist=['x'])
    scorer_clsname = "".join(s.capitalize() for s in name.split('_'))
    return getattr(scorer_module, scorer_clsname)


################################################################################
# Scorer class
################################################################################

class Scorer(object):

    USE_DAT_MATRIX_AND_DISTRIBUTION = False
    USE_SIM_MATRIX = False
    USE_BG_DISTRIBUTION = False
    SKIP_ADJUSTMENTS = False

    def __init__(self, dat_file=None, s_matrix_file=None, bg_distribution_file=None):
        """
        - sim_matrix: the similarity (scoring) matrix to be used. Not all
          methods will use this parameter.
        - bg_distribution: a list containing an amino acid probability distribution. Not
          all methods use this parameter. The default is the blosum62 background, but
          other distributions can be given.
        """
        if self.USE_DAT_MATRIX_AND_DISTRIBUTION and \
                (self.USE_SIM_MATRIX or self.USE_BG_DISTRIBUTION):
            raise Exception()
        # Data file from Kosiol & Goldman 04
        # http://www.ebi.ac.uk/goldman/dayhoff/
        if self.USE_DAT_MATRIX_AND_DISTRIBUTION:
            # TODO: specify dat file in arguments
            self.sim_matrix, self.bg_distribution = read_dat_file(dat_file)
        # Code from Capra & Singh 07
        if self.USE_SIM_MATRIX:
            self.sim_matrix = read_scoring_matrix(s_matrix_file)
        # Code from Capra & Singh 07
        if self.USE_BG_DISTRIBUTION:
            self.bg_distribution = get_distribution_from_file(bg_distribution_file)

    def score(self, alignment, window_size, window_lambda,
            gap_cutoff, gap_penalty, normalize_scores, **kwargs):
        """
        - alignment: Alignment object
        - window_size: set to >1 for window scoring
        - window_lambda: for window method linear combinatioe
        - gap_cutoff: columns with >gap_cutoff gaps are not scored
        - gap_penalty: penalize gaps by this amount
        - normalize_scores: let scores have mean 0, stdev 1 for this alignment
        """
        scores = []
        for i in range(len(alignment.msa[0])):
            col = get_column(i, alignment.msa)

            if len(col) == len(alignment.msa):
                if self.SKIP_ADJUSTMENTS or gap_percentage(col) <= gap_cutoff:
                    scores.append(self.score_col(col, alignment.seq_weights, gap_penalty, alignment))
                else:
                    #XXX: used to be -1000. I adjusted to 0 so the output looks nicer, I think?
                    scores.append(0)
            else:
                sys.stderr.write("Missing sequences in column %d\n" % i)
        if not self.SKIP_ADJUSTMENTS and window_size > 0:
            scores = window_score(scores, window_size, window_lambda)
        if not self.SKIP_ADJUSTMENTS and normalize_scores:
            scores = calc_z_scores(scores, -999)
        return scores

    def score_col(self, col, seq_weights, gap_penalty=1, alignment=None):
        """
        - col: the column to be scored.
        - seq_weights: an array of floats that is used to weight the contribution
          of each seqeuence. If the len(seq_weights) != len(col), then every sequence
          gets a weight of one.
        - gap_penalty: a binary variable: 0 for no gap penalty and 1
          for gap penalty. The default is to use a penalty. The gap penalty used is
          the score times the fraction of non-gap positions in the column.
        """
        raise NotImplementedError()


################################################################################
# Scorer inputs
################################################################################

def read_dat_file(dat_file):
    if not dat_file:
        dat_file = "matrix/jtt-dcmut.dat.txt"
    qij = [[]] # first row is empty
    with open(dat_file) as f:
        for line in f:
            line = line.strip()
            if not line:
                break
            qij.append(map(float,line.split()))
        for line in f:
            line = line.strip()
            if line:
                distr = map(float,line.split())
                break
    for i, row in enumerate(qij):
        row.append(0)
        col = (qij[j][i] for j in xrange(i+1,len(qij)))
        row += col
        row[i] = -sum(row)
    qij = np.matrix(qij)
    distr = np.matrix(distr).T
    return qij, distr


def read_scoring_matrix(sm_file):
    """
    Read in a scoring matrix from a file, e.g., blosum80.bla, and return it
    as an array.
    """
    if not sm_file:
        sm_file = "matrix/blosum62.bla"

    aa_index = 0
    first_line = 1
    row = []
    list_sm = [] # hold the matrix in list form

    try:
        matrix_file = open(sm_file, 'r')
        for line in matrix_file:
            if line[0] != '#' and first_line:
                first_line = 0
                if len(amino_acids) == 0:
                    for c in line.split():
                        aa_to_index[string.lower(c)] = aa_index
                        amino_acids.append(string.lower(c))
                        aa_index += 1

            elif line[0] != '#' and first_line == 0:
                if len(line) > 1:
                    row = line.split()
                    list_sm.append(row)
        matrix_file.close()
    except IOError, e:
        print "Could not load similarity matrix: %s. Using identity matrix..." % sm_file
        # Return identity matrix
        for i in xrange(20):
            row = []
            for j in xrange(20):
                if i == j:
                    row.append(1)
                else:
                    row.append(0)
            list_sm.append(row)
        return list_sm

    # if matrix is stored in lower tri form, copy to upper
    if len(list_sm[0]) < 20:
        for i in range(0,19):
            for j in range(i+1, 20):
                list_sm[i].append(list_sm[j][i])

    for i in range(len(list_sm)):
        for j in range(len(list_sm[i])):
            list_sm[i][j] = float(list_sm[i][j])

    return list_sm


def get_distribution_from_file(fname):
    """
    Read an amino acid distribution from a file. The probabilities should
    be on a single line separated by whitespace in alphabetical order as in
    amino_acids above. # is the comment character.
    """
    if not fname:
        fname = "matrix/blosum62.distribution"

    distribution = []
    try:
        f = open(fname)
        for line in f:
            if line[0] == '#': continue
            line = line[:-1]
            distribution = line.split()
            distribution = map(float, distribution)
        f.close()

    except IOError, e:
        print e, "Using default (BLOSUM62) background."
        return []

    # use a range to be flexible about round off
    if .997 > sum(distribution) or sum(distribution) > 1.003:
        print "Distribution does not sum to 1. Using default (BLOSUM62) background."
        print sum(distribution)
        return []

    return distribution


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
