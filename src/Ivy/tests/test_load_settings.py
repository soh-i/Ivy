import unittest
from Ivy.load_settings import Setting


class TestSettings(unittest.TestCase):
    def setUp(self):
        self.setting = Setting()

    def test_load(self):
        ivy_conf = self.setting.load('IVY_SETTINGS')
        exp_len = 3
        self.assertEqual(len(ivy_conf.keys()), exp_len)

        edit_bench_conf = self.setting.load('EDIT_BENCH_SETTINGS')
        exp_len = 2
        self.assertEqual(len(edit_bench_conf.keys()), exp_len)
        
    def test_pprint(self):
        self.assertIsNone(self.setting.pprint('EDIT_BENCH_SETTINGS'))


if __name__ == '__main__':
    unittest.main()
    
        
