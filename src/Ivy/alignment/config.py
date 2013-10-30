from Ivy.utils import ImutableDict
from collections import namedtuple

class AlignmentConfig_old(object):
    '''
    AlignmentConfig class provides to generate/access namedtupled filtering params
    >>> cnf = AlignmentConfig()
    >>> print cnf.config.mapq
    >>> 25
    '''
    
    def __init__(self):
        self.conf = self.__init_conf()
        
    def __init_conf(self):
        __conf = namedtuple(
            'Conf', [
                'is_duplicate',
                'is_unmapped',
                'is_deletion',
                'is_proper_pair',
                'is_qcfail',
                'is_secondary',
                'mate_is_reverse',
                'mate_is_unmapped',
                'mapq',
                'base_qual',
                'edit_ratio',
                'edit_base_c',
                'mutation_type_c',
            ])
        __c = __conf(
            is_duplicate=False,
            is_unmapped=False,
            is_deletion=False,
            is_proper_pair=True,
            is_qcfail=False,
            is_secondary=True,
            mate_is_reverse=True,
            mate_is_unmapped=False,
            mapq=25,
            base_qual=25,
            edit_ratio=0.1,
            edit_base_c=10,
            mutation_type_c=1,
        )
        return __c
        
