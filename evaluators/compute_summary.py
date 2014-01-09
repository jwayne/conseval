from __future__ import division
import numpy as np
from evaluate import get_batchscores


def compute_summary(dataset_name):
    n_sites = []
    n_positives = []
    n_seqs = []
    n_seqs_orig = []
    for alignment, _ in get_batchscores(dataset_name):
        n_sites.append( len(alignment.msa[0]) )
        n_positives.append( alignment.testset.count(1) )
        n_seqs.append( len(alignment.msa) )
        n_seqs_orig.append( alignment.orig_num_sequences )

    print "Avg # seqs per alignment: %d" % np.mean(n_seqs)
    print "Avg # seqs per alignment before filtering: %d" % np.mean(n_seqs_orig)
    print "Avg # sites per alignment: %d" % np.mean(n_sites)
    print "Avg %% positives per alignment: %f" % np.mean(np.array(n_positives) / np.array(n_sites))
    print "%% positives total: %f" % ( np.sum(np.array(n_positives)) / np.sum(np.array(n_sites)) )
