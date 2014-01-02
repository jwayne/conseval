import random
import os

from params import Params, ParamDef, WithParams
from phylotree import get_phylotree, read_phylotree, check_phylotree
from seqweights import get_seq_weights
from utils.bio import iupac_alphabet, get_column


class Alignment(WithParams):

    params = Params(
        ParamDef("test_file", None,
            help="path to test file of scores to use for this alignment"),
        ParamDef("parse_testset_fn", None,
            help="function to parse testset fields. Meaningful only if test_file set"),
        ParamDef("tree_file", None,
            help="path to custom phylogenetic tree to use for this alignment"),
    )

    MAX_SEQUENCES = 50

    def __init__(self, align_file, **params):
        """
        Loads/calculates input data for `align_file`.  Sets:
        - self.align_file: Filename this alignment was loaded from
        - self.names: Names of sequences in self.msa
        - self.msa: List of equal-length lists of nucleotides/AAs
        - self.testset: Test labels of each column, if available

        self.msa is cleaned to have no more than self.MAX_SEQUENCES sequences. This
        is so phylogenetic tree calculation does not take too long.

        self.msa is cleaned so that the first sequence contains no gaps, and so that
        no column contains all gaps.
        """
        super(Alignment, self).__init__(**params)

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

        if self.test_file:
            testset = self.parse_testset_fn(self.test_file, msa[0])
        else:
            testset = None

        inds = []
        for i in xrange(len(msa[0])):
            if msa[0][i] == '-':
                continue
            col = get_column(i, msa)
            if col.count('-') == len(col):
                continue
            inds.append(i)
        msa = [[row[i] for i in inds] for row in msa]
        if testset:
            testset = [testset[i] for i in inds]

        self.align_file = os.path.abspath(align_file)
        self.names = names
        self.msa = msa
        self.testset = testset
        self._phylotree = None
        self._seq_weights = None


    def get_phylotree(self, n_bootstrap=0, overwrite=False):
        """
        Phylogenetic tree computed on the alignment.

        Caches the loaded/computed phylogenetic tree after the first call. Setting `overwrite`
        ignores the cached value and re-computes the tree.
        """
        if overwrite and self.tree_file:
            raise ValueError("Cannot overwrite tree for alignment when given tree_file")
        if not self._phylotree or overwrite:
            if self.tree_file:
                tree = read_phylotree(self.tree_file)
                if not check_phylotree(self, tree):
                    raise ValueError("Input tree does not match sequences in alignment")
            else:
                tree = get_phylotree(self, n_bootstrap, overwrite)
            self._phylotree = tree
        return self._phylotree

    def get_seq_weights(self):
        """
        Array of floats that is used to weight the contribution of each
        seqeuence, based on HenikoffHenikoff1994, to emphasize sequences
        that are surprising.  This is only needed when a phylogenetic tree
        is not used; it helps adjust for alignments in which sequences are
        all very similar or very different.
        
        Caches the loaded/computed sequence weights after the first call.
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
