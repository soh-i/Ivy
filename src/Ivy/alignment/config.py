
class AlignmentConfig(object):
    def __init__(self):
        self.conf = self.__set_default()
        
    def __set_default(self):
        initialize = {
            is_duplicate : False,
            is_unmapped : False,
            is_deletion : False,
            is_proper_pair : True,
            is_qcfail : False,
            is_secondary : True,
            mapq : 25,
            mate_is_reverse : True,
            mate_is_unmapped : False,
            base_qual : 25,
        }

        
    def add_conf(self):
        pass

    def read(self):
        pass



