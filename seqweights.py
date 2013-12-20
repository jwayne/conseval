
#####
# Conversion tools
#####

def convert_fname_aln2weights(fname_aln):
    return '.'.join(fname_aln.split('.')[:-1]) + '.weights'


def read_sequence_weights(fname):
    """
    Read in a sequence weight file f and create sequence weight list.
    The weights are in the same order as the sequences each on a new line.
    """
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
    """
    Calculate the sequence weights using the Henikoff '94 method
    for the given msa.
    """

    seq_weights = [0.] * len(msa)
    for i in xrange(len(msa[0])):
        freq_counts = [0] * len(amino_acids)

        col = []
        for j in xrange(len(msa)):
            if msa[j][i] != '-': # ignore gaps
                freq_counts[aa_to_index[msa[j][i]]] += 1

        num_observed_types = 0
        for j in xrange(len(freq_counts)):
            if freq_counts[j] > 0: num_observed_types +=1

        for j in xrange(len(msa)):
            d = freq_counts[aa_to_index[msa[j][i]]] * num_observed_types
            if d > 0:
                seq_weights[j] += 1. / d

    for w in xrange(len(seq_weights)):
        seq_weights[w] /= len(msa[0])

    return seq_weights


def get_seq_weights(fname_aln, msa):
    """
    Get the appropriate sequence weights for `fname_aln`.  If they don't exist,
    compute them but do not cache to disk.
    """
    fname_weights = convert_fname_aln2weights(fname_aln)
    seq_weights = read_sequence_weights(fname_weights)
    if not seq_weights or len(seq_weights) != len(msa):
        seq_weights = compute_sequence_weights(msa)
    return seq_weights


if __name__ == "__main__":
    import sys
    from alignment import Alignment
    print get_seq_weights(sys.argv[1], Alignment(sys.argv[1]).msa)
