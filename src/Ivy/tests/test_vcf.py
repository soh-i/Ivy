try:
    import unittest2 as unittest
except ImportError:
    import unittest
    
import datetime
from Ivy.io.vcf import to_xml, header_column_name, _VCFMetaHeader, VCFWriterDataLine
from Ivy.version import __version__
from Ivy.parse_ivy_opts import CommandLineParser


class TestUtilsFunc(unittest.TestCase):
    def setUp(self):
        pass
        
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
        self.version = __version__
        self.meta_header = _VCFMetaHeader(dict())

    @unittest.skip("CommandLineParser class is not callable on this test")
    def test_build(self):
        self.assertEqual(self.meta_header.build(), "")

    def test__add_date_spec(self):

        d = datetime.datetime.today()
        test_date_info = {'fileformat': 'VCFv4.1',
                          'source': 'Ivy_' + self.version,
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
        self.meta_header = _add_file_spec()

    @unittest.skip('Todo: write test case')
    def test__add_params_spec():
        mock_path = ["file://mock_seq1.fa", "file://mock_seq2.bam"]
        # Todo: write test case
        
class TestVCFHeaderBuilder(unittest.TestCase):
    def setUp(self):
        params = ""
        #builder = VCFHeaderBuilder(params)

    @unittest.skip('Not yet')
    def test_build():
        pass

        
class Test_VCFInfoHeader(unittest.TestCase):
    def setUp(self):
        pass
        
        
class TestVCFWriterDataLine(unittest.TestCase):
    def setUp(self):
        self.mock_data = {'pos': 1102, 'chrom': 'chr12', 'id_': 'ID',
                          'ref': 'A', 'alt': 'G', 'qual': 43, 'filt': False}
        self.mock_vcf = {'coverage': 25, 'dp4': [5, 5, 5, 10], 'allele_freq': 0.92,
                         'ag_freq': 0.38, 'average_baq': 21, 'average_mapq': 43,
                         'strand_bias': 0.002, 'positional_bias': 0.0299,
                         'mismatches': 82, 'matches': 23}
        
        self.w_data_line = VCFWriterDataLine()
    
    def test_write_data_line(self):
        self.assertIsInstance(self.mock_data, dict)
        req_keys = ['pos', 'chrom', 'id_', 'ref', 'alt', 'qual', 'filt']
        for k in req_keys:
            self.assertIn(k, self.mock_data)

    def test__to_info(self):
        _id, val = "test_id", 921
        exp = '{id_:}={val:};'.format(id_=_id, val=val)
        self.assertEqual(self.w_data_line._to_info(_id, val), exp)

    def test_DP(self):
        self.assertIsInstance(self.mock_vcf['coverage'], int)
        self.assertRegexpMatches(self.w_data_line._DP(self.mock_vcf), r'^DP=\d+;$')

    def test_DP4(self):
        self.assertIsInstance(self.mock_vcf['dp4'], list)
        self.assertRegexpMatches(self.w_data_line._DP4(self.mock_vcf), r'^DP4=\d+,\d+,\d+,\d+;$')

    def test_AF(self):
        self.assertIsInstance(self.mock_vcf['allele_freq'], float)
        self.assertRegexpMatches(self.w_data_line._AF(self.mock_vcf), r'^AF=\d+\.\d+;$')

    def test_EF(self):
        self.assertIsInstance(self.mock_vcf['ag_freq'], float)
        self.assertRegexpMatches(self.w_data_line._EF(self.mock_vcf), r'^EF=\d+\.\d+;$')
     
    def test_MAPQ(self):
        self.assertIsInstance(self.mock_vcf['average_mapq'], int)
        self.assertRegexpMatches(self.w_data_line._MAPQ(self.mock_vcf), r'^MAPQ=\d+;$')
     
    def test_BACQ(self):
        self.assertIsInstance(self.mock_vcf['average_baq'], int)
        self.assertRegexpMatches(self.w_data_line._BACQ(self.mock_vcf), r'^BACQ=\d+;$')

    def test_SB(self):
        self.assertIsInstance(self.mock_vcf['strand_bias'], float)
        self.assertRegexpMatches(self.w_data_line._SB(self.mock_vcf), r'^SB=\d+\.\d+;$')
     
    def test_PB(self):
        self.assertIsInstance(self.mock_vcf['positional_bias'], float)
        self.assertRegexpMatches(self.w_data_line._PB(self.mock_vcf), r'^PB=\d+\.\d+;$')

    def test_MIS(self):
        self.assertIsInstance(self.mock_vcf['mismatches'], int)
        self.assertRegexpMatches(self.w_data_line._MIS(self.mock_vcf), r'^MIS=\d+;$')
        
    def test_MA(self):
        self.assertIsInstance(self.mock_vcf['matches'], int)
        self.assertRegexpMatches(self.w_data_line._MA(self.mock_vcf), r'^MA=\d+;$')

    def test_write_info_line(self):
        # raises AssertionError, e.g. DP=25; != coverage=25;
        
        s =  self.w_data_line._DP(self.mock_vcf)
        s += self.w_data_line._DP4(self.mock_vcf)
        s += self.w_data_line._AF(self.mock_vcf)
        s += self.w_data_line._EF(self.mock_vcf)
        s += self.w_data_line._MAPQ(self.mock_vcf)
        s += self.w_data_line._BACQ(self.mock_vcf)
        s += self.w_data_line._SB(self.mock_vcf)
        s += self.w_data_line._PB(self.mock_vcf)
        s += self.w_data_line._MA(self.mock_vcf)
        s += self.w_data_line._MIS(self.mock_vcf)
        
        exp = str()
        for content in self.mock_vcf.items():
            if isinstance(content[1], list):
                # to string
                exp += self.w_data_line._to_info(content[0], ",".join(str(_) for _ in content[1]))
            elif isinstance(content[1], str) or isinstance(content[1], int) or isinstance(content[1], float):
                exp += self.w_data_line._to_info(content[0], content[1])
        self.assertNotEqual(s, exp)
        

if __name__ == '__main__':
    unittest.main(verbosity=2)
    
