import os.path
import collections

def find_app_root():
    '''
    Absolute path to your project root from setup.py location
    '''
    root = os.path.dirname(__file__)
    while not os.path.exists(os.path.join(root, 'setup.py')):
        root = os.path.abspath(os.path.join(root, os.path.pardir))
    return root

def __end_url_basename(p):
    """Returns the final component of a pathname"""
    i = p.rfind('/') + 1
    return p[i:]

    
class ImutableDict(collections.Mapping):
    '''
    ImutableDict class generates imutalbe dictionary
    >>> conf = {"A": 192, "test": 293, "hoo": 93}
    >>> imutable = ImutableDict(conf)
    >>> for _ in imutable:
            print _ , imutable[_],
    >>> A 192 test 293 hoo 93
    >>> imutable['A'] = 'AA'
    >>> TypeError: ImutableDict object does not support item assignment
    >>> im.replace('A', 209) # create NEW imutable dictionary
    >>> print im
    ImutableDict(A=209, test=293, hoo=93)
    '''
    
    def __init__(self, dic):
        self.__dict = dict(dic)
        self.hash = None

    def __getitem__(self, key):
        return self.__dict[key]

    def __len__(self):
        return len(self.__dict)

    def __iter__(self):
        return iter(self.__dict)

    def __hash__(self):
        if self.hash is None:
            self.hash = hash(frozenset(self.__dict.items()))
        return self.hash

    def __eq__(self, other):
        return self.__dict == other.__dict

    def __setitem__(self, key, value):
        raise TypeError, ('%s object does not support item assignment'
                          % (self.__class__.__name__))
        
    def __repr__(self):
        args = [
            '{key}={value}'.format(key=key, value=value)
            for key, value in self.__dict.items()
            ]
        args_str = '(' + ', '.join(args) + ')'
        return self.__class__.__name__ + args_str
        
    def replace(self, key, value):
        self.__dict = ImutableDict.__add__({key:value}).__dict

        
if __name__ == '__main__':
    c = {"A":212, 'B':2192}
    im = ImutableDict(c)
    c = {"A":212, 'B':212}
    print im
