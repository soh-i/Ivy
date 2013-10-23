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

    def make_header(self):
        return self.__info_header()
        
    def __info_header(self):
        prefix = '##INFO='
        info_h = ''
        for k in self.params:
            info_h += ','.join([prefix + '<ID=' + k,
                              'Number=' + '1',
                              'Type='+'Interger',
                              'Description=' + '"' + k + 'filtering params' + '">\n'])
        return info

    def info(self):
        info = ''
        for k in self.params:
              info += ';'.join([k + "=" + str(self.params[k]), ''])
        return info
        
            
    def format(self):
        for i in range(10):
            yield i

    def filter(self):
        pass

    def spec(self):
        return self.__spec
        
    
if __name__ == '__main__':
    writer = VCFWriter()
    print writer.info()
    

