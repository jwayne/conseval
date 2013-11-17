#!/usr/bin/env python

################################################################################
# score_conservation.py - Copyright Tony Capra 2007 - Last Update: 03/09/11
#
# 03/09/11 - default window size = 3, as stated in help
# 06/21/09 - seq specific output now compatible with ConCavity
# 06/21/09 - numarray only included when vn_entropy is used
# 08/15/08 - added property_relative_entropy scoring method
# 08/15/08 - added equal sequence length check
# 01/07/08 - added z-score normalization option (-n)
# 01/07/08 - added seq. specific output option (-a)
# 11/30/07 - read_scoring_matrix now returns list rather than array.
# 11/30/07 - added window lambda command line option (-b)
# 07/05/07 - fixed gap penalty cutoff (<=) and error message
# 06/26/07 - fixed read_clustal_align bug
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA
#
# -----------------------------------------------------------------------------
#
# Dependencies:
# numarray (for von Neumann entropy method) -
# http://sourceforge.net/project/showfiles.php?group_id=1369&release_id=223264
#
#
# See usage() for usage.
#
# This program supports the paper: Capra JA and Singh M.
# Predicting functionally important residues from sequence
# conservation. Bioinformatics. 23(15): 1875-1882, 2007.
# Please cite this paper if you use this code.
#
# It contains code for each method of scoring conservation that is
# evaluated: sum of pairs, weighted sum of pairs, Shannon entropy,
# Shannon entropy with property groupings (Mirny and Shakhnovich 95,
# Valdar and Thornton 01), relative entropy with property groupings
# (Williamson 95), von Neumann entropy (Caffrey et al 04), relative
# entropy (Samudrala and Wang 06), and Jensen-Shannon divergence
# (Capra and Singh 07).
#
# The code distributed with Mayrose et al 04 is used for Rate4Site. As
# of today it can be obtained from:
# http://www.tau.ac.il/~itaymay/cp/rate4site.html
#
#
# See models.scorer.Scorer for the function signature for the scoring
# functions.
#
################################################################################

import math, os, sys, getopt
from utils import *
from models.scorer import Scorer


def read_fasta_alignment(filename):
    """ Read in the alignment stored in the FASTA file, filename. Return two
    lists: the identifiers and sequences. """

    names = []
    alignment = []
    cur_seq = ''

    f = open(filename)
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
            alignment.append(cur_seq.replace('B', 'D').replace('Z', 'Q').replace('X', '-'))
            cur_seq = ''
        elif line[0] in iupac_alphabet:
            cur_seq += line.replace('\r', '')
    f.close()

    # add the last sequence
    cur_seq = cur_seq.upper()
    for i, aa in enumerate(cur_seq):
        if aa not in iupac_alphabet:
            cur_seq = cur_seq.replace(aa, '-')
    alignment.append(cur_seq.replace('B', 'D').replace('Z', 'Q').replace('X', '-'))

    return names, alignment


def read_clustal_alignment(filename):
    """ Read in the alignment stored in the CLUSTAL file, filename. Return
    two lists: the names and sequences. """

    names = []
    alignment = []

    f = open(filename)
    for line in f:
        line = line[:-1]
        if len(line) == 0: continue
        if '*' in line: continue
        if 'CLUSTAL' in line: continue

        t = line.split()
        if len(t) == 2 and t[1][0] in iupac_alphabet:
            if t[0] not in names:
                names.append(t[0])
                alignment.append(t[1].upper().replace('B', 'D').replace('Z', 'Q').replace('X', '-').replace('\r', ''))
            else:
                alignment[names.index(t[0])] += t[1].upper().replace('B', 'D').replace('Z', 'Q').replace('X','-').replace('\r', '')
    f.close()

    return names, alignment




################################################################################
# Execution process
################################################################################

DEFAULT_SCORER = "js_divergence"
def get_scorer_cls(name):
    scorer_module = __import__('models.'+name, fromlist=['x'])
    scorer_clsname = "".join(s.capitalize() for s in name.split('_'))
    return getattr(scorer_module, scorer_clsname)

def get_scorer(name, **kwargs):
    try:
        scorer_cls = get_scorer_cls(name)
    except (ImportError, AttributeError), e:
        sys.stderr.write("%s: %s is not a valid scoring method. Using default '%s'\n"
            % (e, arg, DEFAULT_SCORER))
        scorer_cls = get_scorer_cls(DEFAULT_SCORER)
    scorer = scorer_cls(**kwargs)
    return scorer


def run_scorers(align_file, scorers,
        use_seq_weights = True,
        window_size = 3, # 0 = no window
        win_lam = .5, # for window method linear combination
        background_name = 'blosum62',
        gap_cutoff = .3,
        gap_penalty = 1,
        normalize_scores = False,
        ):

    # Ingest inputs - alignment
    try:
        names, alignment = read_clustal_alignment(align_file)
        if not names:
            names, alignment = read_fasta_alignment(align_file)
    except IOError, e:
        raise IOError("%s. Could not find %s. Exiting..." % (e, align_file))

    if len(alignment) != len(names) or alignment == []:
        raise ValueError("Unable to parse alignment.")
    seq_len = len(alignment[0])
    for i, seq in enumerate(alignment):
        if len(seq) != seq_len:
            raise ValueError("Sequences of different lengths: %s (%d) != %s (%d)."
                % (names[0], seq_len, names[i], len(seq)))

    # Ingest inputs - weights?
    seq_weights = []
    if use_seq_weights:
        align_suffix = align_file.split('.')[-1]
        seq_weights = load_sequence_weights(align_file.replace('.%s' % align_suffix, '.weights'))
        if not seq_weights:
            seq_weights = calculate_sequence_weights(alignment)
    if len(seq_weights) != len(alignment):
        # Unclear if you should raise an error if seq weights sucks and use_seq_weights is True...
        seq_weights = [1.] * len(alignment)

    # calculate scores
    if isinstance(scorers, Scorer):
        # for backwards compat
        scorers = [scorers]
    all_scores = [] #list of scores by each scorer
    for scorer in scorers:
        scores = scorer.score(alignment,
                window_size, win_lam,
                seq_weights,
                gap_cutoff, gap_penalty,
                normalize_scores)
        all_scores.append(scores)

    return zip(*all_scores), alignment, names


def print_scores(scores, alignment, names,
        scorer_names,
        out_fname = None,
        seq_specific_output = None, # name of sequence if True
        use_seq_weights = True,
        window_size = 3, # 0 = no window
        win_lam = .5, # for window method linear combination
        background_name = 'blosum62',
        gap_penalty = 1,
        normalize_scores = False,
        ):

    #####
    # TODO: print for each scorer
    #####
    if out_fname:
        try:
            f = open(out_fname, 'w')
        except IOError, e:
            sys.stderr.write("Could not open %s for output. Printing results to standard out.\n" % out_fname)
    else:
        f = sys.stdout

    f.write("# %s -- %s - window_size: %d - window lambda: %.2f - "
            "background: %s - seq. weighting: %s - gap penalty: %d - normalized: %s\n"
            % (align_file, type(scorer).__name__, window_size, win_lam,
            background_name, use_seq_weights, gap_penalty, normalize_scores))

    # handle print of output relative to specific sequence
    ref_seq_num = None
    if seq_specific_output:
        if seq_specific_output not in names:
            sys.stderr.write("Sequence %s not found in alignment. Using default output format...\n" % seq_specific_output)
            seq_specific_output = 0
        else:
            ref_seq_num = names.index(seq_specific_output)

    if seq_specific_output:
        f.write("# reference sequence: %s\n" % seq_specific_output)
        f.write("# align_column_number\tamino acid\n%s\n" % "\t".join(scorer_names))
    else:
        f.write("# align_column_number\tcolumn\n%s\n" % "\t".join(scorer_names))

    # print scores
    scores = zip(scores)
    for i, score in enumerate(scores):
        if ref_seq_num is not None:
            cur_aa = get_column(i, alignment)[ref_seq_num]
            if cur_aa == '-': continue
            annotation = cur_aa
        else:
            annotation = "".join(get_column(i, alignment))
        f.write("%d\t%s\t%s\n" % (i, annotation, "\t".join(("%.5f" % s for s in score))))

    if out_fname:
        f.close()



################################################################################
# Cmd line driver
################################################################################

USAGE = """\nUSAGE:\npython score_conservation.py [options] alignfile\n\t -alignfile must be in fasta or clustal format."""

def parse_args():
    # TODO: usage thing
    parser = argparse.ArgumentParser(description=USAGE)
    parser.add_argument('align_file')
    parser.add_argument('-o', dest='out_fname', type=str, default=None,
        help="name of output file. Default=output to screen [filename]")
    parser.add_argument('-l', dest='use_seq_weights', type='store_false',
        help="use sequence weighting. Default=True [True|False]")
    parser.add_argument('-p', dest='gap_penalty', type=int, default=1,
        help="gap penalty. Lower the score of columns that contain gaps. Default=1")
    parser.add_argument('-m', dest='s_matrix_file', type=str, default=None,
        help="similarity matrix file, *.bla or *.qij. Default=matrix/blosum62.bla [filename]")
    parser.add_argument('-d', dest='bg_distribution_file', type=str, default=None,
        help="background distribution file, e.g., swissprot.distribution. Default=BLOSUM62 background [filename]")
    parser.add_argument('-w', dest='window_size', type=int, default=3,
        help="window size. Number of residues on either side included in the window. Default=3 [int]")
    # TODO: check 0 <= win_lam <= 1
    parser.add_argument('-b', dest='window_lambda', type=float, default=.5,
        help="lambda for window heuristic linear combination. Default=.5 [real in [0,1]]")
    # TODO: check 0 <= win_lam <= 1
    parser.add_argument('-g', dest='gap_cutoff', type=float, default=.3,
        help="gap cutoff. Do not score columns that contain more than gap cutoff fraction gaps. Default=.3 [real in [0, 1)]")
    parser.add_argument('-a', dest='seq_specific_output', type='store_true',
        help="reference sequence. Print scores in reference to a specific sequence (ignoring gaps). Default prints the entire column. [sequence name]")
    parser.add_argument('-n', dest='normalize_scores', type='store_true',
        help="normalize scores. Print the z-score (over the alignment) of each column raw score. Default=False")
    # TODO: get append right
    parser.add_argument('-s', dest='scorer_names', type='store_append',
        help="conservation estimation method")

    args = parser.parse_args()
    return args


def main(args=None):
    align_file = args.align_file

    # Can we convert this to a dict easily?
    kwargs = args

    scorer_names = args.scorer_names
    if not scorer_names:
        scorer_names.append(DEFAULT_SCORER)
    scorers = [get_scorer(name, **kwargs) for name in scorer_names]

    scores, alignment, names= run_scorers(align_file, scorers, **kwargs)

    print_scores(scores, alignment, names, **kwargs)


if __name__ == "__main__":
    main()
