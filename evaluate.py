#!/usr/bin/python
import argparse
import os
import sys


def get_evaluator(name):
    try:
        ev_module = __import__('evaluators.'+name, fromlist=['x'])
        fn_name = name.split('.')[-1]
        fn = getattr(ev_module, fn_name)
    except (ImportError, AttributeError), e:
        raise ImportError("%s: %s is not a valid evaluator." % (e, name))
    return fn


################################################################################
# Cmd line driver
################################################################################

def main():
    parser = argparse.ArgumentParser(
        description="Run evaluators of scoring methods using the output from batch running the scorers on a dataset(s)")

    parser.add_argument('evaluator_name',
        help="name of evaluator")
    parser.add_argument('batch_output_dirs',
        help="directory(s) containing scores, i.e. output directory(s) from batch running scorers")

    args = parser.parse_args()

    name = args.evaluator_name
    evaluator = get_evaluator(name)
    args.scores_dir


if __name__ == "__main__":
    main()
