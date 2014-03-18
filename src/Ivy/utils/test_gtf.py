try:
    import unittest2 as unittest
except ImportError:
    import unittest

from Ivy.utils.gtf import GTF
from Ivy.base import Utils
import os.path
import pysam


class TestGTF(unittest.TestCase):
    def setUp(self):
        self.path_to_sample = os.path.join(Utils.find_app_root(), 'sample')
        self.path_to_mock_data = os.path.join(Utils.find_app_root(), 'sample/genes.gtf')
        self.gtf = GTF(self.path_to_mock_data)
        
    def test__prepare(self):
        self.assertTrue(os.path.isdir(self.path_to_sample))
        self.assertTrue(os.path.isfile(self.path_to_mock_data + '.gz'))

    def test_fetch_gtf(self):
        for data in self.gtf.fetch_gtf(contig="chr2L", end=9839):
            self.assertIsInstance(data, pysam.TabProxies.GTFProxy)
            self.assertRegexpMatches(str(data.contig), r'^.+$')
            self.assertRegexpMatches(str(data.start), r'^\d+$')
            self.assertRegexpMatches(str(data.end), r'^\d+$')
            self.assertRegexpMatches(data.strand, r'^[+-]$')
            
    def test_gene_in_region(self):
        for gene in self.gtf.gene_in_region(contig="chr2L", end=9839):
            self.assertIsInstance(gene, str)
            self.assertRegexpMatches(gene, r'^[\d\w]*.*')

    def test_subset_of_feature_in_region(self):
        for sub in self.gtf.subset_of_feature_in_region(contig="chr2L", end=9839,
                                                        types=["transcript_id", 'tss_id']):
            self.assertEqual(len(sub), 2)
            self.assertIsInstance(sub, dict)
            self.assertRegexpMatches(sub.keys()[1], r'^transcript_id$')
            self.assertRegexpMatches(sub.keys()[0], r'^tss_id$')
            

    def test_strand_info(self):
        self.assertRegexpMatches(self.gtf.strand_info(contig="chr2L", start=9839), r'^[-]$')
        for coord in range(2000, 2500):
            self.assertRegexpMatches(self.gtf.strand_info(contig="chr2L", start=coord), r'^[.+-]$')
        

if __name__ == '__main__':
    unittest.main(verbosity=2)
