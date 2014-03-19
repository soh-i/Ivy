import unittest
from Ivy.benchmark.benchmark import *


class TestDarnedDataGenerator(unittest.TestCase):
    def setUp(self):
        self.gen = DarnedDataGenerator(species='human_hg19')

    def test_constractor(self):
        exp = ['human_hg40', 'mice_mm19', 'fly_dm34']
        for _ in exp:
            self.assertRaises(ValueError, DarnedDataGenerator, species=_)

    def test_variable_name(self):
        exp1 = 'human_hg19.txt'
        self.assertEqual(self.gen.filename, exp1)

        exp2 = 'human_hg19'
        self.assertEqual(self.gen.species, exp2)
        
        exp3 = 'http://darned.ucc.ie/static/downloads/hg19.txt'
        self.assertEqual(self.gen.url, exp3)

    def test_darned_to_csv(self):
        pass
 
#    def test_fetch_darned(self):
#        # darned_data_generator = DarnedDataGenerator(species)
#        # self.assertEqual(expected, darned_data_generator.fetch_darned())
#        assert False # TODO: implement your test here
# 
#class TestDarnedReader(unittest.TestCase):
#    def test_db_name(self):
#        # darned_reader = DarnedReader(sp, source)
#        # self.assertEqual(expected, darned_reader.db_name())
#        assert False # TODO: implement your test here
# 
#    def test_path(self):
#        # darned_reader = DarnedReader(sp, source)
#        # self.assertEqual(expected, darned_reader.path())
#        assert False # TODO: implement your test here
# 
#    def test_size(self):
#        # darned_reader = DarnedReader(sp, source)
#        # self.assertEqual(expected, darned_reader.size())
#        assert False # TODO: implement your test here
# 
#    def test_sp(self):
#        # darned_reader = DarnedReader(sp, source)
#        # self.assertEqual(expected, darned_reader.sp())
#        assert False # TODO: implement your test here
# 
#class TestVCFReader(unittest.TestCase):
#    def test_ag_count(self):
#        # v_cf_reader = VCFReader(filename)
#        # self.assertEqual(expected, v_cf_reader.ag_count())
#        assert False # TODO: implement your test here
# 
#    def test_editing_types(self):
#        # v_cf_reader = VCFReader(filename)
#        # self.assertEqual(expected, v_cf_reader.editing_types())
#        assert False # TODO: implement your test here
# 
#    def test_other_mutations_count(self):
#        # v_cf_reader = VCFReader(filename)
#        # self.assertEqual(expected, v_cf_reader.other_mutations_count())
#        assert False # TODO: implement your test here
# 
#    def test_size(self):
#        # v_cf_reader = VCFReader(filename)
#        # self.assertEqual(expected, v_cf_reader.size())
#        assert False # TODO: implement your test here
# 
#    def test_target_count(self):
#        # v_cf_reader = VCFReader(filename)
#        # self.assertEqual(expected, v_cf_reader.target_count(types))
#        assert False # TODO: implement your test here
# 
#    def test_vcf_name(self):
#        # v_cf_reader = VCFReader(filename)
#        # self.assertEqual(expected, v_cf_reader.vcf_name())
#        assert False # TODO: implement your test here
# 
#class test__CSVReader(unittest.TestCase):
#    def test_ag_count(self):
#        # ___csv_reader = __CSVReader(filename)
#        # self.assertEqual(expected, ___csv_reader.ag_count())
#        assert False # TODO: implement your test here
# 
#    def test_name(self):
#        # ___csv_reader = __CSVReader(filename)
#        # self.assertEqual(expected, ___csv_reader.name())
#        assert False # TODO: implement your test here
# 
#    def test_other_mutations_count(self):
#        # ___csv_reader = __CSVReader(filename)
#        # self.assertEqual(expected, ___csv_reader.other_mutations_count())
#        assert False # TODO: implement your test here
# 
#    def test_size(self):
#        # ___csv_reader = __CSVReader(filename)
#        # self.assertEqual(expected, ___csv_reader.size())
#        assert False # TODO: implement your test here
# 
#class TestBenchmark(unittest.TestCase):
#    def test_f_measure(self):
#        # benchmark = Benchmark(answer, predict)
#        # self.assertEqual(expected, benchmark.f_measure())
#        assert False # TODO: implement your test here
# 
#    def test_precision(self):
#        # benchmark = Benchmark(answer, predict)
#        # self.assertEqual(expected, benchmark.precision())
#        assert False # TODO: implement your test here
# 
#    def test_recall(self):
#        # benchmark = Benchmark(answer, predict)
#        # self.assertEqual(expected, benchmark.recall())
#        assert False # TODO: implement your test here

if __name__ == '__main__':
    unittest.main(verbosity=2)
