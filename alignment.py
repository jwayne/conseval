import random
from phylotree import get_tree
from utils import aa_to_index, amino_acids, iupac_alphabet


class Alignment(object):

    MAX_SEQUENCES = 50

    def __init__(self, align_file):
        """
        Loads/calculates input data for `align_file`.  Sets:
        - self.msa: List of equal-length lists of nucleotides/AAs
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

        # Filter alignment if too many sequences
        if len(msa) > self.MAX_SEQUENCES:
            import ipdb
            ipdb.set_trace()
            random.seed(1000)
            inds = random.sample(range(1,len(msa)), self.MAX_SEQUENCES-1)
            names = [names[0]] + [names[ind] for ind in inds]
            msa = [msa[0]] + [msa[ind] for ind in inds]

        self.align_file = align_file
        self.names = names
        self.msa = msa
        self._phylotree = None
        self._seq_weights = None

    def get_phylotree(self):
        """
        Cached phylogenetic tree.
        """
        if not self._phylotree:
            self._phylotree = get_tree(self.align_file)
        return self._phylotree

    def get_seq_weights(self):
        """
        Cached sequence weights for each column.
        an array of floats that is used to weight the contribution of each
        seqeuence. If the len(seq_weights) != len(col), then every sequence
        gets a weight of one.
        """
        if not self._seq_weights:
            self._seq_weights = get_seq_weights(self.align_file, self.msa)
        return self._seq_weights


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
