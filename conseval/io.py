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


def list_scorer_params(*scorers):
    out = []
    out.append("# Timestamp: %s" % get_timestamp())
    for scorer in scorers:
        out.append("# Scorer: %s" % scorer.name)
        params = scorer.get_params()
        for k, v in params:
            out.append("# \t%s: %s" % (k, v))
    return "\n".join(out) + "\n"
