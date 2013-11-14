import os.path
import collections
import Ivy.version

def die(msg=''):
    raise SystemExit(msg)

class Utils(object):
    '''
    Utility class, for example, 
    >>> Utils.find_app_root()
    >>> /Users/yukke/dev/Ivy
    '''
    
    @staticmethod
    def find_app_root():
        '''Absolute path to your project root from setup.py location'''
        root = os.path.dirname(__file__)
        while not os.path.exists(os.path.join(root, 'setup.py')):
            root = os.path.abspath(os.path.join(root, os.path.pardir))
        return root

    @staticmethod
    def end_url_basename(p):
        """Returns the final component of a pathname"""
        i = p.rfind('/') + 1
        return p[i:]


class AttrDict(dict):
    '''
    How to use a dot . to access members of diction" in StackOverFlow
    http://stackoverflow.com/a/12187277
    '''
    
    def __init__(self, d=None, create=True):
        if d is None:
            d = {}
        supr = super(AttrDict, self)
        supr.__setattr__('_data', d)
        supr.__setattr__('__create', create)

    def __repr__(self):
        return "%s=>(%s)" % (self.__class__.__name__, dict.__repr__(self._data))

    def __getattr__(self, name):
        try:
            value = self._data[name]
        except KeyError:
            if not super(AttrDict, self).__getattribute__('__create'):
                raise
            value = {}
            self._data[name] = value
                
        if hasattr(value, 'items'):
            create = super(AttrDict, self).__getattribute__('__create')
            return AttrDict(value, create)
        return value
                
    def __setattr__(self, name, value):
        self._data[name] = value
        
    def __getitem__(self, key):
        try:
            value = self._data[key]
        except KeyError:
            if not super(AttrDict, self).__getattribute__('__create'):
                raise
            value = {}
            self._data[key] = value
                
        if hasattr(value, 'items'):
            create = super(AttrDict, self).__getattribute__('__create')
            return AttrDict(value, create)
        return value
                    
    def __setitem__(self, key, value):
        self._data[key] = value
        
    def __iadd__(self, other):
        if self._data:
            raise TypeError("A Nested dict will only be replaced if it's empty")
        else:
            return other
                                                    
if __name__ == '__main__':
    dic = {"fasta": "reference.fasta", "region": { "start": 993, "end": 9199} }
    aa = AttrDict(dic)
    print aa
    print aa.fasta
    aa.region.start = 7423480

    
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
    >>> im.replace('A', 209) # update imutable dictionary
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
        for old_key in self.__dict:
            if old_key == key:
                # replace new value into the dic
                self.__dict[key] = value
        else:
            # Add new element into the original dic
            self.__dict = dict(self.__dict.items() + {key:value}.items())
        return self.__dict
