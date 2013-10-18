import ConfigParser
import os
import sys

class DataConfig(object):
    def __init__(self):
        pass

    def find_app_root(self):
        root = os.path.dirname(__file__)
        while not os.path.exists(os.path.join(root, 'setup.py')):
            root = os.path.abspath(os.path.join(root, os.path.pardir))
        return root

        
        
