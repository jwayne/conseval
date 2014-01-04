#!/usr/bin/python
from alignment import Alignment
from utils.general import get_all_module_names
from scorers.base import get_scorer_cls
import sys

def list_params():
    print "================"
    print "Alignment params:"
    print "================"
    print "\n".join("    %s" % pd for pd in Alignment.params.param_defs)
    print ""
    scorer_names = get_all_module_names('scorers')
    for scorer_name in scorer_names:
        if scorer_name == "base":
            continue
        print "================"
        print "Scorer params for %s:" % scorer_name
        print "================"
        scorer_cls = get_scorer_cls(scorer_name)
        print "\n".join("    %s" % pd for pd in scorer_cls.params.param_defs)
        print ""
    sys.exit(0)

if __name__ == "__main__":
    list_params()
