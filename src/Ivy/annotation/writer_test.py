try:
    import unittest2 as unittest
except ImportError:
    import unittest


class TestUtilsFunc(unittest.TestCase):
    def setUp(self):
        pass
        
    def test_to_xml(self):
        from Ivy.annotation.writer import to_xml
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
        from Ivy.annotation.writer import header_column_name
        self.assertEqual(header_column_name(), '#CHROM\tPOS\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\n')


class Test_VCFMetaHeader(unittest.TestCase):
    def setUp(self):
        from Ivy.annotation.writer import _VCFMetaHeader
        from Ivy.version import __version__
        from Ivy.commandline.parse_ivy_opts import CommandLineParser
        
        self.version = __version__
        self.meta_header = _VCFMetaHeader(dict())

    @unittest.skip("CommandLineParser class is not callable on this test")
    def test_build(self):
        self.assertEqual(self.meta_header.build(), "")

    def test__add_date_spec(self):
        import datetime
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
        from Ivy.annotation.writer import VCFWriterDataLine
        self.mock_data = {'pos': 1102, 'chrom': 'chr12', 'id_': 'ID',
                          'ref': 'A', 'alt': 'G', 'qual': 43, 'filt': False}
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

    @unittest.skip("Not yet")
    def test_write_info_line(self):
        pass

    @unittest.skip("Not yet")
    def test_each_columns(self):
        pass
        


        
if __name__ == '__main__':
    unittest.main(verbosity=9)


    
