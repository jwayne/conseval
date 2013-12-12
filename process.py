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
# See scorer.Scorer for the function signature for the scoring
# functions.
#
################################################################################

import argparse, math, os, sys
from alignment import Alignment
from scorer import get_scorer
from utils import get_column


DEFAULT_SCORER = "js_divergence"


################################################################################
# Execution process
################################################################################

#defaults
def run_scorers(alignment, scorers,
        window_size = 3, # 0 = no window
        window_lambda = .5, # for window method linear combination
        gap_cutoff = .3,
        gap_penalty = 1,
        normalize_scores = False,
        ):
    """
    Returns list of scores for each column, i.e.
        [(score1, score2, score3, ..) for col1,
         (score1, score2, score3, ..) for col2,
         ...
        ]
    """

    # calculate scores
    all_scores = [] #list of scores by each scorer
    for scorer in scorers:
        try:
            scores = scorer.score(alignment,
                    window_size=window_size,
                    window_lambda=window_lambda,
                    gap_cutoff=gap_cutoff,
                    gap_penalty=gap_penalty,
                    normalize_scores=normalize_scores)
        except Exception, e:
            import traceback
            sys.stderr.write("Error scoring %s via %s\n" %
                (alignment.align_file, type(scorer).__name__))
            traceback.print_exc()
            scores = [0] * len(alignment.msa)
        all_scores.append(scores)

    return zip(*all_scores)


#defaults
def print_scores(scores, alignment, scorer_names, params,
        out_fname = None,
        seq_specific_output = None, # name of sequence if True
        ):

    if out_fname:
        try:
            f = open(out_fname, 'w')
        except IOError, e:
            sys.stderr.write("Could not open %s for output. Printing results to standard out.\n" % out_fname)
    else:
        f = sys.stdout

    # Print params
    for pkey in sorted(params.keys()):
        f.write("# %s - %s\n" % (pkey, params[pkey]))

    # handle print of output relative to specific sequence
    ref_seq_num = None
    if seq_specific_output:
        if seq_specific_output not in alignment.names:
            sys.stderr.write("Sequence %s not found in alignment. Using default output format...\n" % seq_specific_output)
            seq_specific_output = 0
        else:
            ref_seq_num = alignment.names.index(seq_specific_output)

    if seq_specific_output:
        f.write("# reference sequence: %s\n" % seq_specific_output)
        f.write("# align_column_number\tamino acid\n%s\n" % "\t".join(scorer_names))
    else:
        f.write("# align_column_number\tcolumn\n%s\n" % "\t".join(scorer_names))

    # print scores
    scores = zip(scores)
    for i, score in enumerate(scores):
        if ref_seq_num is not None:
            cur_aa = get_column(i, alignment.msa)[ref_seq_num]
            if cur_aa == '-': continue
            annotation = cur_aa
        else:
            annotation = "".join(get_column(i, alignment.msa))
        f.write("%d\t%s\t%s\n" % (i, annotation, "\t".join(("%.5f" % s for s in score))))

    if out_fname:
        f.close()



################################################################################
# Cmd line driver
################################################################################

def parse_args():
    #defaults
    parser = argparse.ArgumentParser(description="Produce conservation scores for an alignment file.")
    parser.add_argument('align_file')
    parser.add_argument('-o', dest='out_fname', type=str, default=None,
        help="name of output file. Default=output to screen [filename]")
    parser.add_argument('-l', dest='use_seq_weights', action='store_false',
        help="use sequence weighting. Default=True [True|False]")
    parser.add_argument('-p', dest='gap_penalty', type=int, default=1,
        help="gap penalty. Lower the score of columns that contain gaps. Default=1")
    parser.add_argument('-m', dest='s_matrix_file', type=str, default=None,
        help="similarity matrix file, *.bla or *.qij. Default=matrix/blosum62.bla [filename]")
    parser.add_argument('-d', dest='bg_distribution_file', type=str, default=None,
        help="background distribution file, e.g., swissprot.distribution. Default=BLOSUM62 background [filename]")
    parser.add_argument('-w', dest='window_size', type=int, default=3,
        help="window size. Number of residues on either side included in the window. Default=3 [int]")
    def check_bounds(s):
        f = float(s)
        if 0 <= f <= 1: return f
        else: parser.error("must have 0 <= input (%s) <= 1" % s)
    parser.add_argument('-b', dest='window_lambda', type=check_bounds, default=.5,
        help="lambda for window heuristic linear combination. Default=.5 [real in [0,1]]")
    parser.add_argument('-g', dest='gap_cutoff', type=check_bounds, default=.3,
        help="gap cutoff. Do not score columns that contain more than gap cutoff fraction gaps. Default=.3 [real in [0, 1)]")
    parser.add_argument('-a', dest='seq_specific_output', action='store_true',
        help="reference sequence. Print scores in reference to a specific sequence (ignoring gaps). Default prints the entire column. [sequence name]")
    parser.add_argument('-n', dest='normalize_scores', action='store_true',
        help="normalize scores. Print the z-score (over the alignment) of each column raw score. Default=False")
    parser.add_argument('-s', dest='scorer_names', action='append', default=[],
        help="conservation estimation method")

    args = parser.parse_args()
    return args


def main():
    args = parse_args()
    align_file = args.align_file
    kwargs = vars(args)

    # Get scorer(s)
    scorer_names = args.scorer_names
    scorers = []
    for name in scorer_names:
        scorer = get_scorer(name, args.s_matrix_file, args.bg_distribution_file)
        if scorer:
            scorers.append(scorer)
    if not scorers:
        scorer_names = [DEFAULT_SCORER]
        scorers.append(get_scorer(DEFAULT_SCORER, args.s_matrix_file, args.bg_distribution_file))

    # Get alignment and supplementary material
    alignment = Alignment(align_file, args.use_seq_weights)

    scores = run_scorers(alignment, scorers,
        args.window_size,
        args.window_lambda,
        args.gap_cutoff,
        args.gap_penalty,
        args.normalize_scores)

    print_scores(scores, alignment, scorer_names, kwargs,
        args.out_fname, args.seq_specific_output)


if __name__ == "__main__":
    main()
