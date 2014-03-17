try:
    import unittest2 as unittest
except ImportError:
    import unittest
    
from Ivy.annotation.writer import _VCFMetaHeader
from Ivy.annotation.writer import *
from Ivy.version import __version__


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


class Test_VCFMetaHeader(unittest.TestCase):
    def setUp(self):
        ivy_params = {}
        self.meta_header = _VCFMetaHeader(ivy_params)
        
    @unittest.skip('Not yet')
    def test_build(self):
        self.assertEqual(build(), "")

    #@unittest.skip('Not yet')
    def test__add_date_spec(self):
        import datetime
        d = datetime.datetime.today()
        test_date_info = {'fileformat': 'VCFv4.1',
                          'source': 'Ivy_' + __version__,
                          'filedate': '{0}/{1}/{2} {3:0>2d}:{4:0>2d}:{5:0>2d}'.format(
                              d.year, d.month, d.day, d.hour, d.minute, d.second)
        }
        spec = str()
        for _ in test_date_info:
            spec += "##"+  "=".join(['{meta:s}'.format(
                meta=_), test_date_info[_]]) + "\n"
            
        self.assertIsInstance(test_date_info, dict)
        self.assertEqual(self.meta_header._add_date_spec(), spec)
        self.assertNotEqual(self.meta_header._add_date_spec(), test_date_info)

    @unittest.skip('Not yet')
    def test__add_file_spec():
        pass

    @unittest.skip('Not yet')
    def test__add_params_spec():
        pass

        
class TestVCFHeaderBuilder(unittest.TestCase):
    @unittest.skip('Not yet')
    def setUp(self):
        params = ""
        builder = VCFHeaderBuilder(params)

    @unittest.skip('Not yet')
    def test_build():
        pass

        
class Test_VCFInfoHeader(unittest.TestCase):
    def setUp(self):
        pass
        

        
class TestVCFWriterDataLine(unittest.TestCase):
    def setUp(self):
        pass

        
if __name__ == '__main__':
    unittest.main(verbosity=9)


    
