import os
import random


_DATA_HOME_DIR = os.path.abspath("../../data/cs08/")

class DatasetConfig(object):

    def __init__(self, aln_dir, test_dir, align_to_test, parse_testset_line):
        """
        @param aln_dir:
            directory containing the aln files for this dataset, relative to
            `DATA_HOME_DIR`
        @param test_dir:
            directory containing the test files for this dataset, relative to
            `DATA_HOME_DIR`
        @param align_to_test:
            function to convert the aln filename to the test filename. Directories
            are taken care of automatically
        @param parse_testset_line:
            function to parse testset fields for this dataset
        """
        self.aln_dir = os.path.join(_DATA_HOME_DIR, aln_dir)
        self.test_dir = os.path.join(_DATA_HOME_DIR, test_dir)
        self._align_to_test = align_to_test
        self.parse_testset_line = parse_testset_line

    def get_align_files(self, limit=0):
        """
        Get all aln files for `dataset_name`.
        """
        align_files = []
        for root, _, files in os.walk(self.aln_dir):
            for file in files:
                # TODO: handle extensions besides .aln
                if file.endswith('.aln'):
                    align_file = os.path.join(root, file)
                    # Check if test file exists, if not then don't score this alignment.
                    test_file = self.get_test_file(align_file)
                    if not os.path.exists(test_file):
                        continue
                    align_files.append(align_file)
        if limit and len(align_files) > limit:
            random.seed(1000)
            align_files = random.sample(align_files, limit)
        return align_files

    def get_test_file(self, align_file):
        return os.path.join(self.test_dir,
                self._align_to_test(align_file[len(self.aln_dir)+1:]))


################################################################################
# Dataset configs
################################################################################

def _parse_testset_line_csa(line):
    #val is 1 if catalytic site, -1 if not
    fields = line.split()
    val = int(int(fields[1]) > 0)
    pos = int(fields[0])
    return pos, val

def _parse_testset_line_ec(line):
    #val is distance of site from ligand atom
    fields = line.split()
    val = float(fields[2])
    pos = int(fields[0])
    aa = fields[1]
    return pos, val, aa

DATASET_CONFIGS = {
    'csa': DatasetConfig(
        aln_dir = 'conservation_alignments/csa_hssp',
        test_dir = 'cat_sites',
        align_to_test = lambda x: x[:-3] + 'cat_sites',
        parse_testset_line = _parse_testset_line_csa,
    ),
    'ec': DatasetConfig(
        aln_dir = 'conservation_alignments/ec_hssp',
        test_dir = 'lig_distance',
        align_to_test = lambda x: os.path.join(os.path.dirname(x), os.path.split(x)[-1][:6]+ '.dist_to_lig'),
        parse_testset_line = _parse_testset_line_ec,
    ),
}
