from __future__ import division
import numpy as np
import os

from conseval.params import ParamDef
from conseval.utils.bio import amino_acids


N_STATES = 20
PRECISION = .1**5


class SubstitutionModel(object):
    """
    Class to store an amino acid substitution model of an instantaneous
    rate matrix Q and stationary distribution PI, and to calculate the
    probability matrix P(t) = e^(Qt) for any time t>0.

    Code by Josh Chen 2013.
    """

    def __init__(self, dat_file):
        S = [[]] # first row is empty, as the file doesn't contain diagonals
        with open(dat_file) as f:
            for line in f:
                line = line.strip()
                if not line:
                    if len(S) < N_STATES:
                        continue
                    else:
                        break
                row = map(float,line.split())
                if len(row) != len(S):
                    raise ValueError("Bad input rate matrix: Must be in PAML format")
                S.append(row)
            for line in f:
                line = line.strip()
                if line:
                    self.freqs = np.array(map(float,line.split()))
                    if len(self.freqs) != N_STATES:
                        raise ValueError("Bad input equilibrium distribution: Must be in PAML format")
                    break
        # Fill out upper triangle so matrix S is symmetric.
        for j, row in enumerate(S):
            row.append(0)
            col = (S[i][j] for i in xrange(j+1,len(S)))
            row += col
            # Fill in diagonals so that row sums of Q are 0, i.e.
            # sum_j q_ij = sum_j s_ij * pi_j = 0
            row[j] = - np.dot(row, self.freqs) / self.freqs[j]

        self.S = np.matrix(S)
        self.PI = np.diag(self.freqs)

        # Calculate Q = S * PI according to comment in dat file at
        # http://www.ebi.ac.uk/goldman/dayhoff/jtt-dcmut.dat
        self.Q = self.S * self.PI
        # Q should be normalized so that sum_i sum_(j!=i) q_ij * freqs_i = 1
        # according to KosiolGoldman2005 pg 1 + the same comment above.  Verify this.
        Q_rowsums = np.asarray(np.sum(self.Q,1)).reshape(-1) - np.diag(self.Q)
        Q_rowprod = np.dot(self.freqs, Q_rowsums)
        if np.abs(Q_rowprod - 1) > PRECISION:
            raise ValueError("sum_i sum_(j!=i) q_ij * freqs_i = %f != 1" % Q_rowprod)

        # Calculate A = PI^(1/2) * S * PI^(1/2) according to SchabauerEtAl2012 pg 4
        PI_pow_poshalf = np.diag(np.diag(self.PI)**.5)
        PI_pow_neghalf = np.diag(np.diag(PI_pow_poshalf)**(-1))
        self.A = PI_pow_poshalf * self.S * PI_pow_poshalf
        # Calculate Eigvecs * diag(Eigvals) * Eigvecs.T = A according to
        # MolerVan_Loan2003 pg 20, SchabauerEtAl2012 pg 4.
        self.A_eigvals, A_eigvecs = np.linalg.eig(self.A)
        if np.any(np.abs(A_eigvecs * np.diag(self.A_eigvals) * A_eigvecs.T - self.A) > PRECISION):
            raise ValueError("A_eigvecs * diag(A_eigvals) * A_eigvecs.T != A")
        # Prepare multipliers for calc_P.  See calc_P docs for details.
        self.calc_P_left = PI_pow_neghalf * A_eigvecs
        self.calc_P_right = A_eigvecs.T * PI_pow_poshalf
        # Verify that probability calculation is OK.
        P = self.calc_P()
        if np.any(np.abs(self.freqs*P - self.freqs) > PRECISION):
            raise ValueError("pi * P != pi")

    def calc_P(self, t=1):
        """
        Compute the probability matrix P
            P = PI^(-1/2) * e^(At) * PI^(1/2)
            [src: MolerVan_Loan2003 pg 20, SchabauerEtAl2012 pg 4]
        where
            e^(At) = Eigvecs * e^(Eigvals*t) * Eigvecs^-1
            [src: SchabauerEtAl2012 pg 4]
        for Eigvecs, Eigvals of A.
        Note that to speed computation, we pre-compute
            self.calc_P_left = PI^(-1/2) * Eigvecs
            self.calc_P_right = Eigvecs^-1 * PI^(1/2)

        The P returned contains elements P_ij = P(this=j|parent=i)
        """
        A_eigvals_exp = np.diag(np.exp(self.A_eigvals*t))
        return self.calc_P_left * A_eigvals_exp * self.calc_P_right


def read_sim_matrix(sm_file):
    """
    Read in a scoring matrix from a file, e.g., blosum80.bla, and return it
    as an array.
    
    Code by Tony Capra 2007.
    """
    aa_index = 0
    first_line = 1
    row = []
    list_sm = [] # hold the matrix in list form
    with open(sm_file) as matrix_file:
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
                    row = map(float, line.split())
                    list_sm.append(row)
    # if matrix is stored in lower tri form, copy to upper
    if len(list_sm[0]) < 20:
        for i in range(0,19):
            for j in range(i+1, 20):
                list_sm[i].append(list_sm[j][i])
    return list_sm


def read_bg_distribution(fname):
    """
    Read an amino acid distribution from a file. The probabilities should
    be on a single line separated by whitespace in alphabetical order as in
    amino_acids above. # is the comment character.

    Code by Tony Capra 2007.
    """
    distribution = None
    with open(fname) as f:
        for line in f:
            if line[0] == '#': continue
            if distribution:
                raise ValueError("Distribution file has multiple lines")
            fields = line.strip().split()
            distribution = np.array(map(float, fields))
    # use a range to be flexible about round off
    if .997 > sum(distribution) or sum(distribution) > 1.003:
        raise ValueError("Distribution sums to %f != 1." % sum(distribution))
    # sanity check no amino acid has 0 frequency
    if np.any(distribution == 0):
        raise ValueError("Distriution gives some amino acid(s) 0 frequency")
    return distribution


paramdef_sub_model = ParamDef(
        'sub_model', 'sub_models/lg_LG.PAML.txt', os.path.abspath,
        load_fxn=SubstitutionModel,
        help=".dat matrix file of rate matrix and bg distribution")
paramdef_sim_matrix = ParamDef(
        'sim_matrix', 'sub_models/blosum62.bla', os.path.abspath,
        load_fxn=read_sim_matrix,
        help="similarity matrix file, *.bla or *.qij")
paramdef_bg_distribution = ParamDef(
        'bg_distribution', 'sub_models/blosum62.distribution', os.path.abspath,
        load_fxn=read_bg_distribution,
        help="background distribution file, e.g., swissprot.distribution")

def _helper_optional(fn):
    def _clean(x):
        if x: return fn(x)
        return None
    return _clean
paramdef_bg_distribution_optional = ParamDef(
        'bg_distribution', 'sub_models/blosum62.distribution', _helper_optional(os.path.abspath),
        load_fxn=_helper_optional(read_bg_distribution),
        help="background distribution file, e.g., swissprot.distribution")
