"""
INTREPID method of scoring (Sankararaman and Sjolander 08).
Computes JS divergence, for the column, for the subtree rooted at each
node in the path from the root of a phylogenetic tree to the leaf node
of interest (i.e. the target protein), and returns:

   score = max_S ( JS(S,x) - avg_x(S,x) )

where S is a subtree, and x is the column.

Code by Josh Chen 2013.
"""
from __future__ import division
import numpy as np

from scorers.base import Scorer
from scorers.cs07.js_divergence import JsDivergence
from substitution import paramdef_bg_distribution
from utils.bio import get_column


class Intrepid(Scorer):

    params = Scorer.params.with_defaults({
        'normalize': True
    }).extend(paramdef_bg_distribution)

    def __init__(self, **params):
        super(Intrepid, self).__init__(**params)

        params = {
            "bg_distribution": self._bg_distribution,
            "gap_cutoff": 1,
            "use_gap_penalty": False,
            "use_seq_weights": False,
            "window_size": 0,
            "normalize": False,
        }
        self.subscorer = JsDivergence(**params)


    def _score(self, alignment):
        tree = alignment.get_phylotree()
        names_map = dict((name, i) for i, name in enumerate(alignment.names))

        # Find path in tree from root to first sequence in alignment
        # Note that path does not include the root node
        p_name = alignment.names[0]
        p_node = None
        terminals = tree.get_terminals()
        for node in terminals:
            if node.name == p_name:
                p_node = node
                break
        assert p_node
        path = tree.get_path(p_node)

        # Get sequences in each subtree in the path.
        msas = [alignment.msa]
        for subtree in path:
            msas.append([alignment.msa[names_map[node.name]] for node in subtree.get_terminals()])
        seq_weightss = [[1] * len(msa) for msa in msas]

        # subscores = [[score for each subtree] for each col]
        subscores = np.zeros((len(alignment.msa[0]), len(msas)))
        for subtree_ind in xrange(len(msas)):
            msa = msas[subtree_ind]
            seq_weights = seq_weightss[subtree_ind]
            for col_ind in xrange(len(msa[0])):
                col = get_column(col_ind, msa)
                subscores[col_ind][subtree_ind] = self.subscorer._score_col(col, seq_weights)
        avg_subscores = np.mean(subscores,0)
        subscores -= avg_subscores

        scores = list(np.max(subscores,1))
        return scores
