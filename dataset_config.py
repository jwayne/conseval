import os
import random


_DATA_HOME_DIR = os.path.abspath("../../data/cs08/")

class DatasetConfig(object):

    def __init__(self, aln_dir, test_dir, align_to_test, parse_testset_fn):
        """
        @param aln_dir:
            directory containing the aln files for this dataset
        @param test_dir:
            directory containing the test files for this dataset
        @param align_to_test:
            function to convert the aln filename to the test filename
        @param parse_testset_fn:
            function to parse testset fields for this dataset
        """
        self.aln_dir = os.path.join(_DATA_HOME_DIR, aln_dir)
        self.test_dir = os.path.join(_DATA_HOME_DIR, test_dir)
        self.align_to_test = align_to_test
        self.parse_testset_fn = parse_testset_fn

    def get_filenames(self, limit=0):
        """
        Get all aln files for `dataset_name`.
        """
        filenames = []
        for root, _, files in os.walk(self.aln_dir):
            for file in files:
                #TODO: not just aln
                if file.endswith('.aln'):
                    # Get the 'tail' of the path after aln_dir
                    align_file = os.path.join(root, file)
                    align_tail = align_file[len(self.aln_dir)+1:]
                    # Check if test file exists, if not then don't score this alignment.
                    test_file = os.path.join(self.test_dir, self.align_to_test(align_tail))
                    if not os.path.exists(test_file):
                        continue
                    # Add (align_file - aln_dir) to list of filenames
                    filenames.append(align_tail)
        if limit and len(filenames) > limit:
            random.seed(1000)
            filenames = random.sample(filenames, limit)
        return filenames


################################################################################
# Dataset configs
################################################################################

def _parse_testset_csa_fields(fields):
    #result is 1 if catalytic site, 0 if not
    result = int(int(fields[1]) > 0)
    pos = int(fields[0])
    return pos, result

def _parse_testset_ec_fields(fields):
    #result is distance of site from ligand atom
    result = float(fields[2])
    pos = int(fields[0])
    return pos, result

DATASET_CONFIGS = {
    'csa': DatasetConfig(
        aln_dir = 'conservation_alignments/csa_hssp',
        test_dir = 'cat_sites',
        align_to_test = lambda x: x[:-3] + 'cat_sites',
        parse_testset_fn = _parse_testset_csa_fields,
    ),
    'ec': DatasetConfig(
        aln_dir = 'conservation_alignments/ec_hssp',
        test_dir = 'lig_distance',
        align_to_test = lambda x: os.path.join(os.path.dirname(x), os.path.split(x)[-1][:6]+ '.dist_to_lig'),
        parse_testset_fn = _parse_testset_ec_fields,
    ),
}
