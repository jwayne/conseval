#!/usr/bin/python
"""
Create an Alignment object for each alignment file in a dataset.
This can be used to easily test parsing of a test set.
"""
from alignment import Alignment
from datasets import DATASET_CONFIGS as DC
import sys

dc=DC[sys.argv[1]]
alns=[]
align_files=dc.get_align_files()

print len(align_files)
for align_file in align_files:
    test_file = dc.get_test_file(align_file)
    aln = Alignment(align_file, test_file=test_file, parse_testset_fn=dc.parse_testset_fn)
