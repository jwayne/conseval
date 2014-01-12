import matplotlib.pyplot as plt
import numpy as np
import random
import os
from conseval.alignment import Alignment, MockAlignment
from conseval.datasets import DATASET_CONFIGS
from evaluate import get_batchscores, get_batchscore_dir
from conseval.io import read_batchscores
from conseval.utils.bio import get_column, amino_acids
from conseval.utils.general import weighted_choice
from scorers.rate4site_eb import Rate4siteEb, precompute_tree_probs
from scorers.js_divergence import JsDivergence


N_RUNS = 50
def r4s_ppc(dataset_name, **jsd_params):
    afs = list(get_batchscores(dataset_name, align_files_only=True))

    dc = DATASET_CONFIGS[dataset_name]
    r4s_name = 'R4S_EB-vanilla'
    r4s_dir = os.path.join(get_batchscore_dir(dataset_name), r4s_name)

    r4s = Rate4siteEb()
    jsd = JsDivergence(**jsd_params)

    # Choose random alignment/scores pair
    align_file = random.choice(afs)
    test_file = dc.get_test_file(align_file)
    r4s_file = dc.get_out_file(align_file, r4s_dir)
    alignment = Alignment(align_file, test_file=test_file, parse_testset_fn=dc.parse_testset_fn)
    n_seqs = len(alignment.msa)
    n_sites = len(alignment.msa[0])

    fig = plt.figure()
    ax = plt.gca()
    inds = range(n_sites)

    rates = read_batchscores(r4s_file)
    tree = alignment.get_phylotree()
    root = tree.root
    names_map = dict((name,i) for i,name in enumerate(alignment.names))

    # Pre-compute probabilities for branch, for every site (i.e., rate)
    P_cached = precompute_tree_probs(tree, rates, r4s.sub_model)
    for r in P_cached:
        for node in P_cached[r]:
            # We treat root separately
            if node != root:
                cum_probs = np.cumsum(P_cached[r][node],axis=1)
                if not np.all(np.abs(cum_probs[:,-1]-1) < 1e-4):
                    raise ValueError("Bad probability matrix")
                cum_probs[:,-1] = 1
                P_cached[r][node] = np.array(cum_probs)
    root_freqs = np.cumsum(r4s.sub_model.freqs)
    if not abs(root_freqs[-1]-1) < 1e-4:
        raise ValueError("Bad probability matrix")
    root_freqs[-1] = 1

    # Repeat replication for N_RUNS
    jsd_rep_scores_all = []
    for n in xrange(N_RUNS):
        # For each site, generate amino acids for each sequence using that site's rate
        # DFS through tree, setting amino acids at each node
        msa = [[] for i in xrange(n_seqs)]
        for i in xrange(n_sites):
            r = rates[i]
            aa_ind = weighted_choice(root_freqs)
            bfs = [(aa_ind,node) for node in root.clades]
            for aa_ind, node in bfs:
                aa_ind = weighted_choice(P_cached[r][node][aa_ind])
                if node.is_terminal():
                    msa[names_map[node.name]].append(amino_acids[aa_ind])
                else:
                    bfs += ((aa_ind,child) for child in node.clades)
        aln_rep = MockAlignment(alignment.names, msa, tree, alignment.get_seq_weights)
        jsd_rep_scores = jsd.score(aln_rep)
        jsd_rep_scores_all.append(jsd_rep_scores)
        ax.scatter(inds, jsd_rep_scores, color='k', alpha=0.2)

    jsd_orig_scores = jsd.score(alignment)
    ts = alignment.testset
    pos_inds = [i for i in inds if ts[i]]
    pos_orig_scores = [jsd_orig_scores[i] for i in pos_inds]
    neg_inds = [i for i in inds if not ts[i]]
    neg_orig_scores = [jsd_orig_scores[i] for i in neg_inds]
    ax.scatter(pos_inds, pos_orig_scores, color='g')
    ax.scatter(neg_inds, neg_orig_scores, color='r')
    ax.set_ylabel('JSD score')

    rep_scores_per_col = np.array(jsd_rep_scores_all).T
    means = np.mean(rep_scores_per_col, axis=1)
    stds = np.std(rep_scores_per_col, axis=1)
    zscores = (np.array(jsd_orig_scores) - means) / stds
    ax2 = ax.twinx()
    ax2.plot(inds, zscores, '--')
    ax2.set_ylabel('Deviations')

    plt.xlim(inds[0], inds[-1])
    plt.xlabel('Sites')

    plt.figure()
    plt.scatter(rates, zscores)
    plt.xlabel('Rates')
    plt.xlabel('Deviations')

    plt.show(block=False)
    import ipdb
    ipdb.set_trace()
