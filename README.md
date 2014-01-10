# ConsEval

Simplifying the implemention, evaluation, and running of sequence conservation scoring methods for
locating protein functional sites.

# Setup

The following dependencies are required:
* Python 2.6+
* numpy
* BioPython ([biopython.org], or `sudo easy_install -f http://biopython.org/DIST/ biopython`)
* PyYAML for parsing config files ([pyyaml.org], or `sudo easy_install pyyaml`)
* PhyML for computing phylogenetic trees ([code.google.com/p/phyml/])

You should also create the following folders (or symlinks) in the main folder (i.e., the folder with the scripts `score.py`, `batchscore.py`, and `evaluate.py` in it):
* `input/`, for all the input data files (alignments and test labels).  I personally symlink this to an external data folder, whose subfolders are each datasets of alignments and test labels.  I also create config files for `batchscore.py` (see next section) here.
* `output/`, for all the output that will be created by `batchscore.py`.

This code should work on any UNIX-based OS (Linux or MAC).


# Getting started

This code supports three functions.

First, you can compute conservation rates, using a specified scoring method,
on a single alignment stored in CLUSTAL or FASTA format.  Note that available
scoring methods can be found in `scorers/`.

```
# Estimate conservation rates using the JS divergence scoring method.
./score.py cs07.js_divergence examples/1dup_A_hssp-filtered.aln

# Estimate conservation rates using the JS divergence scoring method
# but with different parameters.
./score.py cs07.js_divergence examples/1dup_A_hssp-filtered.aln -p lambda_prior=0.1 -p gap_cutoff=0 -p window_size=0

# Estimate conservation rates using the JS divergence scoring method
# but with a custom phylogenetic tree for the alignment.
./score.py cs07.js_divergence examples/1dup_A_hssp-filtered.aln -a tree_file=examples/1dup_A_hssp-filtered.tree

# List parameters for the JS divergence scoring method.
./score.py cs07.js_divegence -l

# List parameters for all available scoring methods.
./score.py -l

# Help with scoring single alignment files.
./score.py -h
```


Second, you can run a (parallelized) batch job of estimating conservation rates, on multiple
datasets of alignments (paired with test labels) and using multiple scoring
methods.  (Test labels are needed because I assume that you are running a batch
job with the purpose generate scores for evaluating different scoring methods.
This requirement can be removed in the future if one wishes to run batch jobs to
score unlabeled datasets, without the purpose of evaluating different scoring
methods.)

```
# Run a batch scoring job, specified by the YAML config file at examples/example.yaml.
./batchscore.py examples/example.yaml
```


Third, you can run evaluations of the output from your batch jobs.
Available evaluators can be found in `evaluators/`.

```
# Draw PR and ROC curves for the output of batch jobs on the 'csa' dataset
# for scoring runs 
./evaluate pr_roc csa 
```


# Customization

It should be easy to:
* Define and implement new scoring methods in `scorers/`
* Define and implement new evaluation methods in `evaluators/`
* Add different substitution models to `sub_models/`
It is perhaps easiest to look at the other files in these folders as examples
for how such scorers, evaluators, and substitution models should be written/formatted.


# Acknowledgements

This project is the result of part of my 2013 senior fall independent work at
Princeton.  Special thanks to Mona Singh for being my advisor.

Some code was initially based off of Tony Capra's work 2007-2011,
available at [http://compbio.cs.princeton.edu/conservation/].
The datasets this code is intended to be run on can also be found at the above link.
The paper behind that code and the datasets support the following paper:
* Capra JA and Singh M. Predicting functionally important residues from
sequence conservation. Bioinformatics. 23(15):1875-82, 2007.  

The substitution matrices included in this package use files packaged as part
of the code above, as well as files from the following papers:
* Kosiol C and Goldman N. Different Versions of the Dayhoff Rate Matrix.
Molecular Biology and Evolution. 22(2):193-199, 2005.
* Le S.Q., Gascuel O. LG: An Improved, General Amino-Acid Replacement
Matrix. Molecular Biology and Evolution. 25(7):1307-20, 2008.

Happy scoring!
