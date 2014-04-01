import os.path
import pprint
from Ivy.settings import (
    IVY_SETTINGS,
    EDIT_BENCH_SETTINGS
)


class Setting(object):
    def __init__(self):
        self.ivy = IVY_SETTINGS
        self.edit_bench = EDIT_BENCH_SETTINGS
        
    def load(self, cls):
        if cls == 'IVY_SETTINGS':
            return self.ivy
        elif cls == 'EDIT_BENCH_SETTINGS':
            return self.edit_bench

    def pprint(self, cls):
        pp = pprint.PrettyPrinter(indent=1)
        if cls == 'IVY_SETTINGS':
            return pp.pprint(self.ivy)
        elif cls == 'EDIT_BENCH_SETTINGS':
            return pp.pprint(self.edit_bench)
        else:
            raise KeyError("Do NOT match given your key named '{0}'".format(cls))
            
        
        
        
        
if __name__ == '__main__':
    setting = Setting()
    conf =  setting.load('IVY_SETTINGS')
    setting.pprint('EDIT_BENCH_SETTINGS')
