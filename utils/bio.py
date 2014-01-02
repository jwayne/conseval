from __future__ import division
import numpy as np


amino_acids = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V', '-']
iupac_alphabet = ["A", "B", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "U", "V", "W", "Y", "Z", "X", "*", "-"]

# dictionary to map from amino acid to its row/column in a similarity matrix
aa_to_index = dict((aa,i) for i,aa in enumerate(amino_acids))


def get_column(col_num, msa):
    """Return the `col_num`-th column of `msa` as a list."""
    return [seq[col_num] for seq in msa]


################################################################################
# Frequency Count and Gap Penalty
################################################################################

PSEUDOCOUNT = .0000001

def weighted_freq_count_pseudocount(col, seq_weights, pc_amount):
    """
    Return the weighted frequency count for a column--with pseudocount.
    """
    freq_counts = np.zeros(len(amino_acids))
    for j,aa in enumerate(col):
        freq_counts[aa_to_index[aa]] += seq_weights[j]
    freq_counts /= np.sum(freq_counts)
    return freq_counts


def weighted_gap_penalty(col, seq_weights):
    """
    Calculate the simple gap penalty multiplier for the column. If the
    sequences are weighted, the gaps, when penalized, are weighted
    accordingly.
    """
    gap_sum = sum(seq_weights[i] for i in xrange(len(col)) if col[i] == '-')
    return 1 - (gap_sum / sum(seq_weights))
