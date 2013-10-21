from Ivy.benchmark import *
import unittest

class VCFReaderTest(unittest.TestCase):
    
    def setUp(self):
        self.db = VCFReader('../../../data/test_data.vcf')
        
    def test_count(self):
        count = self.db.cnt()
        self.assertEqual(count, 13298)

    def test_name(self):
        name = self.db.vcf_name()
        self.assertEqual(name, 'test_data.vcf')

    def test_ag_count(self):
        ag = self.db.ag_count()
        self.assertEqual(ag, 1609)

    def test_other_mutation_count(self):
        c = self.db.other_mutations_count()
        self.assertEqual(c, 11689)
        
    def suite(self):
        suite = unittest.TestSuite()
        suite.addTest(unittest.makeSuite(VCFReaderTest))
        return suite

if __name__ == '__main__':
    unittest.main(verbosity=2).run(suite)
    
