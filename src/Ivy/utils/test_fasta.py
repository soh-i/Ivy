try:
    import unittest2 as unittest
except ImportError:
    import unittest
from Ivy.utils.fasta import Fasta
from Ivy.base import Utils
import os
import os.path
import pysam


class TestFasta(unittest.TestCase):
    def setUp(self):
        path = 'sample/hg19_mock.fasta'
        self.mock_fasta = os.path.join(Utils.find_app_root(), path)
        self.fasta_inst = Fasta(self.mock_fasta)

    def test_fasta_header(self):
        exp_header = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5']
        unexp_header = [1, 2, 3, 4, 5]
        self.assertEqual(self.fasta_inst.fasta_header(), exp_header)
        self.assertNotEqual(self.fasta_inst.fasta_header(), unexp_header)
                         
    def test_split_by_blocks(self):
        test_num = 2
        self.fasta_inst.split_by_blocks(n=test_num)
        self.assertEqual(len(os.listdir('./block_fasta')), test_num)
        
    def test_chr_size(self):
        exp_len = 5
        self.assertEqual(self.fasta_inst.chr_size(), exp_len)

    def test_split_by_length(self):
        pass

    def test_generate_chrom_blockcs(self):
        exp_size = 3
        result = self.fasta_inst.generate_chrom_blocks(exp_size)
        self.assertEqual(len(result), exp_size)
        self.assertEqual(len(result[0]), 2)
        self.assertEqual(len(result[1]), 2)
        self.assertEqual(len(result[2]), 1)
        

class TestFastaUtils(unittest.TestCase):
    def setUp(self):
        pass

    @unittest.skip("Not yet")
    def test_decode_chr_name_from_file(self):
        pass

    @unittest.skip("Not yet")
    def test_as_single(self):
        pass

    @unittest.skip("Not yet")
    def test_fetch_seq(self):
        pass

    @unittest.skip("Not yet")
    def test_run(self):
        pass

    @unittest.skip("Not yet")
    def test_get_fa_list(self):
        pass


if __name__ == '__main__':
    unittest.main(verbosity=2)
