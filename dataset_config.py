import os
import random


DATA_HOME_DIR = os.path.abspath("../../data/cs07/")
OUT_HOME_DIR = os.path.abspath("../../out/")

class DatasetConfig(object):

    def __init__(self, aln_dir, test_dir, align_to_test, parse_testset_fn):
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
        @param parse_testset_fn:
            function of (test_file, args) to parse the test file for this dataset
        """
        self.aln_dir = os.path.join(DATA_HOME_DIR, aln_dir)
        self.test_dir = os.path.join(DATA_HOME_DIR, test_dir)
        self._align_to_test = align_to_test
        self.parse_testset_fn = parse_testset_fn

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
        return sorted(align_files)

    def get_test_file(self, align_file):
        return os.path.join(self.test_dir,
                self._align_to_test(align_file[len(self.aln_dir)+1:]))

    def get_out_file(self, align_file, out_dir, ext='.res'):
        out_name = ".".join(align_file[len(self.aln_dir)+1:].replace('/', '___').split('.')[:-1])
        return os.path.join(out_dir, out_name + ext)



################################################################################
# Dataset configs
################################################################################

def _parse_testset_csa(test_file, seq):
    """
    Test files have 2-column rows of (position, value).  Positions range from
    0 to len(seq)-1, with all sites in the seq included as rows (even gaps).
    Values are 1 if catalytic site, and -1 if not.
    """
    n_sites = len(seq)
    vals = []
    with open(test_file) as f:
        for line in f:
            if line.startswith('#'):
                continue
            fields = line.split()
            #val is 1 if catalytic site, -1 if not
            pos, val = int(fields[0]), int(int(fields[1]) > 0)
            if pos != len(vals):
                raise ValueError("%r: pos = %d, need %d" % (test_file, pos, len(vals)))
            if pos >= n_sites:
                raise ValueError("%r: pos = %d > len(seq) = %d" % (test_file, pos, n_sites))
            vals.append(val)
    return vals


def _parse_testset_ec(test_file, seq):
    """
    Test files have 3-column rows of (position, aa, value).  Positions are
    difficult to understand (as far as I can tell).  Values are matched to
    the aa; they are the distance of site from a ligand atom.  Each site in
    seq (including gaps) has a row in the test file.

    This is a best-effort match of the test file's sites to the alignment
    sequence's sites. Any unmatches sites in the alignment seq have 'None'
    values in returned list.
    """
    lines = []
    with open(test_file) as f:
        for line in f:
            if line.startswith('#'):
                continue
            fields = line.split()
            pos, aa, val = int(fields[0]), fields[1], float(fields[2])
            lines.append((pos, aa, val))
    seq_nogaps = ''.join([aa for aa in seq if aa != '-'])

    n_sites = len(seq)
    vals = []
    start_pos = None
    for i, (pos, aa, val) in enumerate(lines):
        # Already have vals for every site, so return.
        if len(vals) >= n_sites:
            return vals

        # Gaps in sequence should be matched with None's in the testset.
        seq_aa = seq[len(vals)]
        while seq_aa == '-':
            vals.append(None)
            if len(vals) >= n_sites:
                return vals
            seq_aa = seq[len(vals)]

        # Try to guess sequence alignment via 'pos'
        if start_pos is None:
            start_pos = pos
        curr_pos = pos - start_pos
        if len(vals) > curr_pos:
            # pos skips down incorrectly, so adjust start_pos so pos doesn't skip
            start_pos -= len(vals) - curr_pos
            curr_pos = pos - start_pos

        if seq_aa != aa:
            # The test file is skipping some position indices. This suggests
            # that the test file is missing scores for these indices, but this
            # actually isn't always the case; position indices can be mislabeled
            # in the test file. We want to find out how many sites the test file
            # is actually missing scores for. We do this by finding the first
            # 5-aa match (adjusting length for gaps and for the total length of
            # the sequences) between the test sequence going forward, and the
            # alignment sequence going forward, at various increments (_curr_pos)
            # along the alignment sequence.
            diff_match = -1
            rg = curr_pos - len(vals)
            assert rg >= 0
            rg = max(20, rg)
            for diff in [0, rg] + range(rg):
                _curr_pos = len(vals) + diff
                for j in xrange(5):
                    if _curr_pos+j >= len(seq) or i+j >= len(lines):
                        break
                    _seq_aa = seq[_curr_pos+j]
                    _test_aa = lines[i+j][1]
                    if _seq_aa == '-':
                        break
                    if _seq_aa != _test_aa:
                        diff_match = -1
                        break
                    else:
                        diff_match = _curr_pos - len(vals)
                if diff_match >= 0:
                    break
            if diff_match < 0:
                # Perhaps the test file has an extra row.
                start_pos += curr_pos - len(vals)
                continue
            vals += [None] * diff_match
            start_pos += diff_match
            seq_aa = seq[len(vals)]

        assert seq_aa == aa
        vals.append(val)

    if len(vals) < len(seq):
        vals += [None] * (len(seq) - len(vals))
    return vals


DATASET_CONFIGS = {
    'csa': DatasetConfig(
        aln_dir = 'conservation_alignments/csa_hssp',
        test_dir = 'cat_sites',
        align_to_test = lambda x: x[:-3] + 'cat_sites',
        parse_testset_fn = _parse_testset_csa,
    ),
    'ec': DatasetConfig(
        aln_dir = 'conservation_alignments/ec_hssp',
        test_dir = 'lig_distance',
        align_to_test = lambda x: os.path.join(os.path.dirname(x), os.path.split(x)[-1][:6]+ '.dist_to_lig'),
        parse_testset_fn = _parse_testset_ec,
    ),
}
