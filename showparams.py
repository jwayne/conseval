#!/usr/bin/python
import sys
from conseval.alignment import Alignment
from conseval.scorer import get_scorer_cls
from conseval.utils.general import get_all_module_names


def list_alignment_paramdefs():
    print "================"
    print "Alignment params:"
    print "================"
    print "\n".join("    %s" % pd for pd in Alignment.params.param_defs)
    print ""


def list_scorer_paramdefs():
    scorer_names = get_all_module_names('scorers')
    for scorer_name in scorer_names:
        try:
            scorer_cls = get_scorer_cls(scorer_name)
        except ImportError:
            continue
        print "================"
        print "Scorer params for %s:" % scorer_name
        print "================"
        print "\n".join("    %s" % pd for pd in scorer_cls.params.param_defs)
        print ""


if __name__ == "__main__":
    list_alignment_paramdefs()
    list_scorer_paramdefs()
