try:
    import unittest2 as unittest
except ImportError:
    import unittest
from Ivy.utils.fasta import Fasta
import pysam
    

class TestFasta(unittest.TestCase):
    def setUp(self):
        self.mock_fasta = ''
        self.fasta_inst = Fasta(self.mock_fasta)
        
    def test_fasta_header(self):
        self.assertEqual(self.fasta_inst.fasta_header(), "")
            
    def test_split_by_blocks(self):
        pass

    def test_chr_size(self):
        pass
        
    def test_split_by_length(self):
        pass

    def test_generate_chrom_blockcs(self):
        pass

    def test__merge_list(self):
        pass


class TestFastaUtils(unittest.TestCase):
    def setUp(self):
        pass

    def test_decode_chr_name_from_file(self):
        pass

    def test_as_single(self):
        pass

    def test_fetch_seq(self):
        pass

    def test_run(self):
        pass

    def test_get_fa_list(self):
        pass


if __name__ == '__main__':
    unittest.main(verbosity=2)
