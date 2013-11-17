import math
import os
import sys
from utils import gap_percentage, get_column, amino_acids, aa_to_index


################################################################################
# Scorer inputs
################################################################################

def read_scoring_matrix(sm_file):
    """ Read in a scoring matrix from a file, e.g., blosum80.bla, and return it
    as an array. """
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
    #sim_matrix = array(list_sm, type=Float32)
    #return sim_matrix


def get_distribution_from_file(fname):
    """ Read an amino acid distribution from a file. The probabilities should
    be on a single line separated by whitespace in alphabetical order as in
    amino_acids above. # is the comment character."""

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


def load_sequence_weights(fname):
    """Read in a sequence weight file f and create sequence weight list.
    The weights are in the same order as the sequences each on a new line. """
    seq_weights = []
    try:
        f = open(fname)
        for line in f:
            l = line.split()
            if line[0] != '#' and len(l) == 2:
                seq_weights.append(float(l[1]))
        f.close()
    except IOError, e:
        pass
    return seq_weights


def calculate_sequence_weights(msa):
    """ Calculate the sequence weights using the Henikoff '94 method
    for the given msa. """

    seq_weights = [0.] * len(msa)
    for i in range(len(msa[0])):
        freq_counts = [0] * len(amino_acids)

        col = []
        for j in range(len(msa)):
            if msa[j][i] != '-': # ignore gaps
                freq_counts[aa_to_index[msa[j][i]]] += 1

        num_observed_types = 0
        for j in range(len(freq_counts)):
            if freq_counts[j] > 0: num_observed_types +=1

        for j in range(len(msa)):
            d = freq_counts[aa_to_index[msa[j][i]]] * num_observed_types
            if d > 0:
                seq_weights[j] += 1. / d

    for w in range(len(seq_weights)):
        seq_weights[w] /= len(msa[0])

    return seq_weights


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


################################################################################
# Scorer
################################################################################

# BLOSUM62 background distribution
blosum_background_distr = [0.078, 0.051, 0.041, 0.052, 0.024, 0.034, 0.059, 0.083, 0.025, 0.062, 0.092, 0.056, 0.024, 0.044, 0.043, 0.059, 0.055, 0.014, 0.034, 0.072]

current_dir = os.path.dirname(os.path.realpath(__file__))


class Scorer(object):

    USE_SIM_MATRIX = False
    USE_BG_DISTRIBUTION = False

    def __init__(self, s_matrix_file, bg_distribution_file):
        """
        - sim_matrix: the similarity (scoring) matrix to be used. Not all
          methods will use this parameter.
        - bg_distribution: a list containing an amino acid probability distribution. Not
          all methods use this parameter. The default is the blosum62 background, but
          other distributions can be given.
        """
        if self.USE_SIM_MATRIX:
            if not s_matrix_file:
                s_matrix_file = os.path.join(current_dir, "matrix/blosum62.bla")
            self.sim_matrix = read_scoring_matrix(s_matrix_file)
        if self.USE_BG_DISTRIBUTION:
            # XXX: they used to copy the list.  doesn't seem necessary so i removed it
            self.bg_distribution = blosum_background_distr
            if bg_distribution_file:
                d = get_distribution_from_file(bg_distribution_file)
                if d:
                    sys.stderr.write("WARNING: Could not read bg distribution from %s. "
                                     "Using default bg distribution\n")
                else:
                    self.bg_distribution = d

    def score(self, alignment, window_size, window_lambda, seq_weights,
            gap_cutoff, gap_penalty, normalize_scores):
        """
        - alignment: List of equal-length lists of nucleotides/AAs
        - window_size: set to >1 for window scoring
        - window_lambda: for window method linear combinatioe
        - seq_weights: ?
        - gap_cutoff: columns with >gap_cutoff gaps are not scored
        - gap_penalty: penalize gaps by this amount
        - normalize_scores: let scores have mean 0, stdev 1 for this alignment
        """
        scores = []
        for i in range(len(alignment[0])):
            col = get_column(i, alignment)

            if len(col) == len(alignment):
                if gap_percentage(col) <= gap_cutoff:
                    scores.append(self.score_col(col, seq_weights, gap_penalty))
                else:
                    #XXX: used to be -1000
                    scores.append(0)
            else:
                sys.stderr.write("Missing sequences in column %d\n" % i)
        if window_size > 0:
            scores = window_score(scores, window_size, window_lambda)
        if normalize_scores:
            scores = calc_z_scores(scores, -999)
        return scores

    def score_col(self, col, seq_weights, gap_penalty):
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
