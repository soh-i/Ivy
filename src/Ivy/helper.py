import os.path
import os
import pprint
import collections
import logging
import Ivy.version


class Timer(object):
    def __enter__(self):
        self.start = time.time()
        return self

    def __exit__(self, *args):
        self.end = time.time()
        self.interval = self.end - self.start

        
class IvyLogger(object):
    def __init__(self):
        self.log_fmt = '[%(asctime)s] [%(levelname)s] [%(message)s]'
        logging.basicConfig(level=logging.DEBUG, format=self.log_fmt) #, filename=str(os.getpid()) + "_ivy.log")


def die(msg=''):
    raise SystemExit(msg)

def find_app_root():
    '''Absolute path to your project root from setup.py location'''
    root = os.path.dirname(__file__)
    while not os.path.exists(os.path.join(root, 'setup.py')):
        root = os.path.abspath(os.path.join(root, os.path.pardir))
    return root

def end_url_basename(p):
    """Returns the final component of a pathname"""
    i = p.rfind('/') + 1
    return p[i:]

def convert_base(ref=None, alt=None, strand=None):
    if ref is None or alt is None:
        return False
    dna = {'A': 'T',
           'T': 'A',
           'G': 'C',
           'C': 'G',
           }
    if strand == '+':
        return [ref, alt]
    elif strand == '-':
        pair = []
        _r = dna.get(ref)
        _a = dna.get(alt)
        if _r:
            pair.append(_r)
        else:
            pair.append(ref)
        if _a:
            pair.append(_a)
        # expected: 'A,G'
        else:
            pair.append(alt)
        return pair
    elif strand == '.':
        return [ref, alt]
    else:
        raise ValueError("Invalid strand data. '-/+' is expected, but '{0}' is given".format(strand))

        
