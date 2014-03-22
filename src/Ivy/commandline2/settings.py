import os.path
import ConfigParser
from Ivy.base.Utils import find_app_root

class Setting(object):
    def __init__(self):
        APP_ROOT = Utils.find_app_root()
        FILE = 'conf.yml'
        self.path_to_conf = os.path.join(APP_ROOT, FILE)
        
    def load(self):
        '''
        Return:
         ConfigParser object
        '''
        
        if os.path.isfile(self.path_to_conf):
            config = ConfigParser()
            self.config = config.read(self.path_to_conf)
        else:
            raise ValueError("{0} is not found".format(self.path_to_conf))
    
    def get_path(self):
        return self.path_to_conf

    def show_all(self):
        #self.load()
        pass

if __name__ == '__main__':
    setting = Setting()
    setting.load()
    
