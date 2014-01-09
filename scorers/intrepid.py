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

from conseval.scorer import Scorer, get_scorer_cls
from conseval.params import ParamDef
from conseval.substitution import paramdef_bg_distribution, paramdef_sub_model
from conseval.utils.bio import get_column


class Intrepid(Scorer):

    params = Scorer.params.extend(
        #TODO: clean_fxn shouldn't allow bad subscorer names
        ParamDef('subscorer_cls', 'cs07.js_divergence', load_fxn=get_scorer_cls,
            help="subscorer"),
        ParamDef('normalize_subscores', False, bool,
            help="normalize scores of subscorer for each "),

        #TODO: don't define these here.  not sure how else to do it though
        ParamDef('lambda_pw', .5, float, lambda x: 0<=x<=1,
            help="prior weight lambda_pw in the Jensen-Shannon divergence"),
        paramdef_bg_distribution,
        paramdef_sub_model,
    )

    def __init__(self, **params):
        super(Intrepid, self).__init__(**params)

        subscorer_param_names = set(p.name for p in self.subscorer_cls.params.param_defs)
        # yuck
        for k in params.keys():
            if k not in subscorer_param_names:
                del params[k]
        if 'gap_cutoff' in subscorer_param_names:
            params['gap_cutoff'] = 1
        if 'use_gap_penalty' in subscorer_param_names:
            params['use_gap_penalty'] = False
        if 'use_seq_weights' in subscorer_param_names:
            params['use_seq_weights'] = False
        params.update({
            "window_size": 0,
            "normalize": self.normalize_subscores,
        })
        self.subscorer = self.subscorer_cls(**params)


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



class MockAlignment():
    """
    For sending to Intrepid.subscorer._score(alignment)
    """
    def __init__(self, names, msa, tree, get_seq_weights):
        self.names = names
        self.msa = msa
        self.tree = tree
        self.get_seq_weights = get_seq_weights
    def get_phylotree(self):
        return self.tree
