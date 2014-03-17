try:
    import unittest2 as unittest
except ImportError:
    import unittest
    
from Ivy.annotation.writer import *


class TestUtilsFunc(unittest.TestCase):
    def test_to_xml(self):
        xml_class = "VCF_INFO"
        xml_value = {'id': 'NS', 'num': 1, 'type': 'Integer', 'desc': 'Number of Samples With Data'}
        bad_value = [9212, "hoge", 0.192, list()]
        # passd only dict
        self.assertIsInstance(xml_value, dict)
        for _ in bad_value:
            self.assertNotIsInstance(_, dict)
            
        # test to_xml()
        self.assertEqual(to_xml(xml_class, xml_value),
                         '##VCF_INFO=<ID=NS,Number=1,Type=Integer,Description="Number of Samples With Data">\n')
        
    def test_header_column_name(self):
        self.assertEqual(header_column_name(), '#CHROM\tPOS\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\n')

        
class TestVCFHeaderBuilder(unittest.TestCase):
    def setUp(self):
        pass

        
class Test_VCFInfoHeader(unittest.TestCase):
    def setUp(self):
        pass
        
class Test_VCFMetaHeader(unittest.TestCase):
    def setUp(self):
        pass

        
class TestVCFWriterDataLine(unittest.TestCase):
    def setUp(self):
        pass

        
if __name__ == '__main__':
    unittest.main(verbosity=2)


    
