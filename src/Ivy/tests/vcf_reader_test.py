from Ivy.benchmark import *
import unittest

class VCFReaderTest(unittest.TestCase):
    def setUp(self):
        self.vcf_reader = VCFReader('../../../data/test_data.vcf')
        self.vcf.generate_vcf_set()

    def test_count(self):
        pass

    def test_vcf_name(self):
        pass

    def test_editing_types(self):
        pass

    def test_ag_count(self):
        pass

    def test_other_mutations_count(self):
        pass

    def suite():
        suite = unittest.TestSuite()
        suite.addTest(unittest.makeSuite(VCFReaderTest))
        return suite

if __name__ == '__main__':
    unittest.main(verbosity=2).run(suite)
    
