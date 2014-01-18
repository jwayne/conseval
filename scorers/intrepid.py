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

from conseval.alignment import MockAlignment
from conseval.scorer import Scorer, get_scorer_cls
from conseval.params import ParamDef
from conseval.substitution import paramdef_bg_distribution, paramdef_sub_model
from conseval.utils.bio import get_column


class Intrepid(Scorer):

    params = Scorer.params.extend(
        #TODO: clean_fxn shouldn't allow bad subscorer names
        ParamDef('subscorer_cls', 'js_divergence', load_fxn=get_scorer_cls,
            help="subscorer"),
    )

    def __init__(self, **params):
        super(Intrepid, self).__init__(**params)

        self.subscorer = self.subscorer_cls()


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

        # For each subtree in path, compute scores
        subtree_scores = []
        for subtree in reversed(path):
            # Note that the ordering of sequences in the msa gets messed up.
            inds = [names_map[node.name] for node in subtree.get_terminals()]
            names = [alignment.names[i] for i in inds]
            msa = [alignment.msa[i] for i in inds]
            tree = subtree
            def get_seq_weights():
                if alignment._seq_weights:
                    return alignment._seq_weights
                x=alignment.get_seq_weights()
                return [x[i] for i in inds]
            aln = MockAlignment(names, msa, tree, get_seq_weights)
            subtree_scores.append(self.subscorer.score(aln))
        subtree_scores.append(self.subscorer.score(alignment))

        site_scores = np.array(subtree_scores).T
        site_avgscores = np.mean(site_scores,0)

        scores = list(np.max(site_scores - site_avgscores, 1))
        return scores
