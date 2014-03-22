import os.path
from Ivy.base import Utils
from Ivy.commandline2.settings import IVY_SETTINGS, EDIT_BENCH_SETTINGS

class Setting(object):
    def __init__(self):
        pass
        
    def load(self, cls):
        if cls == 'IVY_SETTINGS':
            return IVY_SETTINGS
        elif cls == 'EDIT_BENCH_SETTINGS':
            return EDIT_BENCH_SETTINGS
        
if __name__ == '__main__':
    setting = Setting()
    print setting.load('EDIT_BENCH_SETTINGS')
    
