from phylotree import get_tree
from utils import aa_to_index, amino_acids, iupac_alphabet


class Alignment(object):

    #defaults
    def __init__(self, align_file, use_seq_weights=True):
        """
        Loads/calculates input data for `align_file`.  Sets:
        - self.msa: List of equal-length lists of nucleotides/AAs
        - self.seq_weights: Weight for each column
        """
        # Ingest alignment
        try:
            names, msa = read_clustal_alignment(align_file)
            if not names:
                names, msa = read_fasta_alignment(align_file)
        except IOError, e:
            raise IOError("%s. Could not find %s. Exiting..." % (e, align_file))

        # Sanity check
        if len(msa) != len(names) or msa == []:
            raise ValueError("Unable to parse alignment.")
        seq_len = len(msa[0])
        for i, seq in enumerate(msa):
            if len(seq) != seq_len:
                raise ValueError("Sequences of different lengths: %s (%d) != %s (%d)."
                    % (names[0], seq_len, names[i], len(seq)))
        if len(set(names)) != len(names):
            raise ValueError("Sequences have duplicate names.")

        # Ingest sequence weights
        seq_weights = []
        if use_seq_weights:
            align_suffix = align_file.split('.')[-1]
            seq_weights = read_sequence_weights(align_file.replace('.%s' % align_suffix, '.weights'))
            if not seq_weights:
                seq_weights = calculate_sequence_weights(msa)
        if len(seq_weights) != len(msa):
            # Unclear if you should raise an error if seq weights sucks and use_seq_weights is True...
            seq_weights = [1.] * len(msa)

        self.align_file = align_file
        self.names = names
        self.msa = msa
        self.seq_weights = seq_weights
        self._phylotree = None

    def get_phylotree(self):
        if not self._phylotree:
            self._phylotree = get_tree(self.align_file)
        return self._phylotree


################################################################################
# Ingesting inputs
################################################################################

def read_fasta_alignment(filename):
    """ Read in the alignment stored in the FASTA file, filename. Return two
    lists: the identifiers and sequences. """
    names = []
    msa = []
    cur_seq = ''
    with open(filename) as f:
        for line in f:
            line = line[:-1]
            if len(line) == 0: continue

            if line[0] == ';': continue
            if line[0] == '>':
                names.append(line[1:].replace('\r', ''))

                if cur_seq != '':
                    cur_seq = cur_seq.upper()
                    for i, aa in enumerate(cur_seq):
                        if aa not in iupac_alphabet:
                            cur_seq = cur_seq.replace(aa, '-')
                msa.append(cur_seq.replace('B', 'D').replace('Z', 'Q').replace('X', '-'))
                cur_seq = ''
            elif line[0] in iupac_alphabet:
                cur_seq += line.replace('\r', '')
    # add the last sequence
    cur_seq = cur_seq.upper()
    for i, aa in enumerate(cur_seq):
        if aa not in iupac_alphabet:
            cur_seq = cur_seq.replace(aa, '-')
    msa.append(cur_seq.replace('B', 'D').replace('Z', 'Q').replace('X', '-'))
    return names, msa


def read_clustal_alignment(filename):
    """Read in the alignment stored in the CLUSTAL file, filename. Return
    two lists: the names and sequences."""
    names = []
    msa = []
    with open(filename) as f:
        for line in f:
            line = line[:-1]
            if len(line) == 0: continue
            if '*' in line: continue
            if 'CLUSTAL' in line: continue

            t = line.split()
            if len(t) == 2 and t[1][0] in iupac_alphabet:
                if t[0] not in names:
                    names.append(t[0])
                    msa.append(t[1].upper().replace('B', 'D').replace('Z', 'Q').replace('X', '-').replace('\r', ''))
                else:
                    msa[names.index(t[0])] += t[1].upper().replace('B', 'D').replace('Z', 'Q').replace('X','-').replace('\r', '')
    return names, msa


def read_sequence_weights(fname):
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
