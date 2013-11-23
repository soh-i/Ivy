import unittest

class TestDarnedDataGenerator(unittest.TestCase):
    def test___init__(self):
        # darned_data_generator = DarnedDataGenerator(species)
        assert False # TODO: implement your test here

    def test_darned_to_csv(self):
        # darned_data_generator = DarnedDataGenerator(species)
        # self.assertEqual(expected, darned_data_generator.darned_to_csv())
        assert False # TODO: implement your test here

    def test_fetch_darned(self):
        # darned_data_generator = DarnedDataGenerator(species)
        # self.assertEqual(expected, darned_data_generator.fetch_darned())
        assert False # TODO: implement your test here

class TestDarnedReader(unittest.TestCase):
    def test___init__(self):
        # darned_reader = DarnedReader(sp, source)
        assert False # TODO: implement your test here

    def test___str__(self):
        # darned_reader = DarnedReader(sp, source)
        # self.assertEqual(expected, darned_reader.__str__())
        assert False # TODO: implement your test here

    def test_db_name(self):
        # darned_reader = DarnedReader(sp, source)
        # self.assertEqual(expected, darned_reader.db_name())
        assert False # TODO: implement your test here

    def test_path(self):
        # darned_reader = DarnedReader(sp, source)
        # self.assertEqual(expected, darned_reader.path())
        assert False # TODO: implement your test here

    def test_size(self):
        # darned_reader = DarnedReader(sp, source)
        # self.assertEqual(expected, darned_reader.size())
        assert False # TODO: implement your test here

    def test_sp(self):
        # darned_reader = DarnedReader(sp, source)
        # self.assertEqual(expected, darned_reader.sp())
        assert False # TODO: implement your test here

class TestVCFReader(unittest.TestCase):
    def test___init__(self):
        # v_cf_reader = VCFReader(filename)
        assert False # TODO: implement your test here

    def test_ag_count(self):
        # v_cf_reader = VCFReader(filename)
        # self.assertEqual(expected, v_cf_reader.ag_count())
        assert False # TODO: implement your test here

    def test_editing_types(self):
        # v_cf_reader = VCFReader(filename)
        # self.assertEqual(expected, v_cf_reader.editing_types())
        assert False # TODO: implement your test here

    def test_other_mutations_count(self):
        # v_cf_reader = VCFReader(filename)
        # self.assertEqual(expected, v_cf_reader.other_mutations_count())
        assert False # TODO: implement your test here

    def test_size(self):
        # v_cf_reader = VCFReader(filename)
        # self.assertEqual(expected, v_cf_reader.size())
        assert False # TODO: implement your test here

    def test_target_count(self):
        # v_cf_reader = VCFReader(filename)
        # self.assertEqual(expected, v_cf_reader.target_count(types))
        assert False # TODO: implement your test here

    def test_vcf_name(self):
        # v_cf_reader = VCFReader(filename)
        # self.assertEqual(expected, v_cf_reader.vcf_name())
        assert False # TODO: implement your test here

class test__CSVReader(unittest.TestCase):
    def test___init__(self):
        # ___csv_reader = __CSVReader(filename)
        assert False # TODO: implement your test here

    def test_ag_count(self):
        # ___csv_reader = __CSVReader(filename)
        # self.assertEqual(expected, ___csv_reader.ag_count())
        assert False # TODO: implement your test here

    def test_name(self):
        # ___csv_reader = __CSVReader(filename)
        # self.assertEqual(expected, ___csv_reader.name())
        assert False # TODO: implement your test here

    def test_other_mutations_count(self):
        # ___csv_reader = __CSVReader(filename)
        # self.assertEqual(expected, ___csv_reader.other_mutations_count())
        assert False # TODO: implement your test here

    def test_size(self):
        # ___csv_reader = __CSVReader(filename)
        # self.assertEqual(expected, ___csv_reader.size())
        assert False # TODO: implement your test here

class TestBenchmark(unittest.TestCase):
    def test___init__(self):
        # benchmark = Benchmark(answer, predict)
        assert False # TODO: implement your test here

    def test___str__(self):
        # benchmark = Benchmark(answer, predict)
        # self.assertEqual(expected, benchmark.__str__())
        assert False # TODO: implement your test here

    def test_f_measure(self):
        # benchmark = Benchmark(answer, predict)
        # self.assertEqual(expected, benchmark.f_measure())
        assert False # TODO: implement your test here

    def test_precision(self):
        # benchmark = Benchmark(answer, predict)
        # self.assertEqual(expected, benchmark.precision())
        assert False # TODO: implement your test here

    def test_recall(self):
        # benchmark = Benchmark(answer, predict)
        # self.assertEqual(expected, benchmark.recall())
        assert False # TODO: implement your test here

if __name__ == '__main__':
    unittest.main()
