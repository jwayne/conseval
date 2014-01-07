"""
Implementation of Henikoff and Henikoff 1994 method to weight sequences
according to their similarity.

Code by Tony Capra 2007.
"""
from __future__ import division
import os
from utils.bio import amino_acids, aa_to_index

#####
# Conversion tools
#####

def get_seq_weights(alignment):
    weights_file = '.'.join(alignment.align_file.split('.')[:-1]) + '.weights'
    seq_weights = read_seq_weights(weights_file)
    if not seq_weights or len(seq_weights) != len(alignment.msa):
        seq_weights = _compute_seq_weights(alignment.msa)
    return seq_weights
        

def read_seq_weights(fname):
    """
    Read in a sequence weight file f and create sequence weight list.
    The weights are in the same order as the sequences each on a new line.
    """
    if not os.path.exists(fname):
        return None
    seq_weights = []
    with open(fname) as f:
        for line in f:
            l = line.split()
            if line[0] != '#' and len(l) == 2:
                seq_weights.append(float(l[1]))
    return seq_weights


def _compute_seq_weights(msa):
    """
    Calculate the sequence weights using the Henikoff '94 method
    for the given msa.
    """
    seq_weights = [0.] * len(msa)
    # For each column
    for i in xrange(len(msa[0])):
        # Find the frequency q of amino acids across all sequences
        freq_counts = [0] * len(amino_acids)
        for j in xrange(len(msa)):
            if msa[j][i] != '-':
                freq_counts[aa_to_index[msa[j][i]]] += 1
        # Find the number of nonzero q's, N
        num_observed_types = 0
        for fc in freq_counts:
            if fc:
                num_observed_types +=1
        # Add 1 / (q_{seq} * N).  This seems kind of weird.
        # 1 / q_{seq} favors sequences with rarer amino acids in the column
        # 1 / N favors adjustments from sites with fewer differences in their column
        for j in xrange(len(msa)):
            if msa[j][i] != '-':
                seq_weights[j] += 1. / (freq_counts[aa_to_index[msa[j][i]]] * \
                                        num_observed_types)
    for w in xrange(len(seq_weights)):
        seq_weights[w] /= len(msa[0])
    return seq_weights


#####

if __name__ == "__main__":
    import sys
    from alignment import Alignment
    print get_seq_weights(sys.argv[1], Alignment(sys.argv[1]).msa)
