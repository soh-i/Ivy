import os.path
from Ivy.base import Utils
from Ivy.commandline2.settings import IVY_SETTINGS, EDIT_BENCH_SETTINGS

class Setting(object):
    def __init__(self):
        self.ivy = IVY_SETTINGS
        self.edit_bench = EDIT_BENCH_SETTINGS
        
    def load(self, cls):
        if cls == 'IVY_SETTINGS':
            return self.ivy
        elif cls == 'EDIT_BENCH_SETTINGS':
            return self.edit_bench
        
if __name__ == '__main__':
    setting = Setting()
    conf =  setting.load('IVY_SETTINGS')
    print conf
    
    

    
