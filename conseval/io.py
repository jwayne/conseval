import os
from conseval.utils.general import get_timestamp


def write_score_helper(score):
    if isinstance(score, float):
        return str(round(score,4))
    elif score is None:
        return '-'
    return str(score)


def read_score_helper(field):
    if field == '-':
        return None
    return float(field)


def write_batchscores(fname, scores):
    if os.path.exists(fname):
        raise IOError("batchscores file %s already exists" % fname)
    with open(fname, 'w') as f:
        f.write("\n".join(map(write_score_helper, scores)))


def read_batchscores(fname):
    with open(fname) as f:
        return map(read_score_helper, f.read().strip().split())


def list_scorer_params(*scorers):
    out = []
    out.append("# Timestamp: %s" % get_timestamp())
    for scorer in scorers:
        out.append("# Scorer: %s" % scorer.name)
        params = scorer.get_params()
        for k, v in params:
            out.append("# \t%s: %s" % (k, v))
    return "\n".join(out) + "\n"
