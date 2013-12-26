import random
import os

from params import Params, ParamDef
from phylotree import get_phylotree, read_phylotree
from seqweights import get_seq_weights
from utils.bio import iupac_alphabet, get_column


class Alignment(object):

    PARAMS = Params(
        ParamDef("test_file", None,
            help="path to test file of scores to use for this alignment"),
        ParamDef("parse_testset_fn", None,
            help="function to parse testset fields"),
        ParamDef("tree_file", None,
            help="path to custom phylogenetic tree to use for this alignment"),
    )

    MAX_SEQUENCES = 50

    def __init__(self, align_file, **params):
        """
        Loads/calculates input data for `align_file`.  Sets:
        - self.msa: List of equal-length lists of nucleotides/AAs
        """
        self.PARAMS.set_params(self, params)

        # Ingest alignment
        try:
            names, msa = read_clustal_alignment(align_file)
            if not names:
                names, msa = read_fasta_alignment(align_file)
        except IOError, e:
            raise IOError("%s. Could not find %s. Exiting..." % (e, align_file))

        # Sanity check
        if not msa:
            raise ValueError("No alignment read.")
        if len(msa) != len(names):
            raise ValueError("Unequal numbers of names and sequences in alignment.")
        seq_len = len(msa[0])
        for i, seq in enumerate(msa):
            if len(seq) != seq_len:
                raise ValueError("Sequences of different lengths: %s (%d) != %s (%d)."
                    % (names[0], seq_len, names[i], len(seq)))
        if len(set(names)) != len(names):
            raise ValueError("Sequences have duplicate names.")

        # Filter alignment if too many sequences.  This is only so that tree
        # computation doesn't take too long.
        self.filtered = False
        self.orig_num_sequences = len(names)
        if self.orig_num_sequences > self.MAX_SEQUENCES:
            self.filtered = True
            random.seed(1000)
            inds = random.sample(range(1,len(msa)), self.MAX_SEQUENCES-1)
            names = [names[0]] + [names[ind] for ind in inds]
            msa = [msa[0]] + [msa[ind] for ind in inds]

        self.align_file = os.path.abspath(align_file)
        self.names = names
        self.msa = msa
        self.testset = None
        self._phylotree = None
        self._seq_weights = None

        if self.test_file:
            self.testset = parse_testset(self.test_file, self.parse_testset_fn, self)

        if self.tree_file:
            # self._phylotree must match msa
            tree = read_phylotree(self.tree_file)
            tree_terminals = tree.get_terminals()
            if len(tree_terminals) != len(self.names) or \
                    set(clade.name for clade in tree_terminals) != set(self.names):
                raise ValueError("Input tree does not match sequences in alignment")
            self._phylotree = tree

    def get_phylotree(self, n_bootstrap=0, overwrite=False):
        """
        Cached phylogenetic tree.
        """
        if not self._phylotree:
            # self._phylotree must match msa
            self._phylotree = get_phylotree(self, n_bootstrap, overwrite)
        return self._phylotree

    def get_seq_weights(self):
        """
        Cached sequence weights for each column.
        an array of floats that is used to weight the contribution of each
        seqeuence. If the len(seq_weights) != len(col), then every sequence
        gets a weight of one.
        """
        if not self._seq_weights:
            self._seq_weights = get_seq_weights(self)
        return self._seq_weights


################################################################################
# Ingesting inputs
################################################################################

def read_fasta_alignment(filename):
    """
    Read in the alignment stored in the FASTA file, filename. Return two
    lists: the identifiers and sequences.
    """
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
    """
    Read in the alignment stored in the CLUSTAL file, filename. Return
    two lists: the names and sequences.
    """
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


def parse_testset(test_file, parse_fields_fn, alignment):
    actual = []
    start_pos = None

    with open(test_file, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            fields = line.split()
            ind, result = parse_fields_fn(fields)
            if start_pos is None:
                start_pos = ind
            pos = ind - start_pos

            if len(actual) > pos:
                # if indices skip down, then re-adjust so it's normal
                start_pos -= len(actual) - pos
            elif pos > len(actual):
                for i in xrange(len(actual), pos):
                    if get_column(i, alignment.msa) == 'X':
                        # if shitty, fill it in
                        actual.append(None)
                        sys.stderr.write("%d: %s\n" % (ind, test_file))
                    else:
                        # otherwise, re-adjust
                        start_pos += 1
            actual.append(result)
    return actual
