import ConfigParser
import os
import sys

def find_app_root():
    root = os.path.dirname(__file__)
    while not os.path.exists(os.path.join(root, 'setup.py')):
        root = os.path.abspath(os.path.join(root, os.path.pardir))
    return root

if __name__ == '__main__':
    print find_app_root()
