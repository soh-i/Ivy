from Ivy.benchmark import *
import unittest

class BenchmarkTest(unittest.TestCase):
    def setUp(self):
        pass

    def tearDown(self):
        pass

    @unittest.skip('not yet')
    def test_sp(self):
        pass

    def suite(self):
        suite = unittest.TestSuite()
        suite.addTest(unittest.makeSuite(BenchmarkTest))
        return suite

if __name__ == '__main__':
    unittest.main(verbosity=2).run(suite)
        
    
