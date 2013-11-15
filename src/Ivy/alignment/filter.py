from Ivy.alignment.stream import AlignmentStream
from Ivy.annotation.writer import VCFWriteHeader

__program__ = 'filter'
__author__ = 'Soh Ishiguro <yukke@g-language.org>'
__license__ = ''
__status__ = 'development'


class AlignmentReadFilter(AlignmentStream):
    def __init__(self):
        super(self.__class__, self).__init__(self)
        #AlignmentStream.__init__(self, self.config.)

class BasicReadFilter(object):
    def __init__(self):
        pass

class StatReadFilter(object):
    def __init__(self):
        pass

class ExtReadFilter(object):
    def __init__(self):
        pass


if __name__ == '__main__':
    ar_filter = AlignmentReadFilter()


    
