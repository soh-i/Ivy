import collections

from alignment import AlignmentConfig


class VCFWriter(object):
    def __init__(self):
        #self.stream = stream
        self.version = '4.1'
        self.__spec = ('http://www.1000genomes.org/wiki/Analysis/'
                       'Variant%20Call%20Format/vcf-variant-call-format-version-41')
        
        p = AlignmentConfig()
        self.params = p.conf

    def info(self):
        prefix = '##INFO='
        for k in self.params:
            yield ','.join([prefix + '<ID=' + k,
                            'Number=' + '1',
                            'Type='+'Interger',
                            'Description=' + '"' + k + 'filtering params' + '">' ])

    def format(self):
        pass

    def header(self):
        self.fileformat = 'VCFv4.1'

        self.reference;
        self.contig;


    def filter(self):
        pass

    def spec(self):
        return self.__spec
        
    
if __name__ == '__main__':
    writer = VCFWriter()
    for i in writer.info():
        print i

