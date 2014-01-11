import matplotlib.pyplot as plt
import numpy as np
from evaluate import get_batchscores, get_out_dir


def compute_variance(dataset_name, *batchscore_ids):
    assert len(batchscore_ids) == 1
    var_list = []
    for alignment, scores_tup in get_batchscores(dataset_name, batchscore_ids):
        var_list.append(np.var(zip(*scores_tup)))
    print np.mean(var_list)
    plt.hist(var_list)
    plt.title('Per-sequence variances of rates r')
    plt.xlabel('Variance')
    plt.ylabel('Count')
    plt.show()
