import unittest
from Ivy.cli.edit_bench_cli_opts import parse_bench_opts
import parser

#class TestEditBenchCliOpts(unittest.TestCase):
#    def __init__(self, testname):
#        super(TestEditBenchCliOpts, self).__init__(testname)
# 
#    def test_cli_parser(self):
#        a = parse_bench_opts()
#        print dir(a)
# 

class TestEditBenchCliOpts(unittest.TestCase):
    def setUp(self):
        self.parser = parse_bench_opts()

    def test(self):
        self.assetEqual(self.parser.sp, 'human_hg19')
    
        
if __name__ == '__main__':
    #suite = unittest.TestSuite()
    #suite.addTest(TestEditBenchCliOpts("test_cli_parser"))
    #unittest.TextTestRunner().run(suite)

    unittest.main()
