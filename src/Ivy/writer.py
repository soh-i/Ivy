import collections
from alignment import AlignmentConfig

class VCFWriter(object):
    def __init__(self, stream):
        self.stream = stream
        self.__spec = {'fileformat': 'VCFv4.1',
                       'source': 'Ivy_v0.0.1',
                       'reference': 'test.fasta',
                       'species': 'Homo Sapiens',
                       'samples': 'sample.bam',
                       
                   }
    
    def merge_filtering_param(self):
        p = AlignmentConfig()
        self.__spec.update({k: p.conf[k] for k in p.conf})
        return self.__spec

    def spec_section(self):
        '''
        spec_section() -> str, return VCF header of basic informations
        e.g. ##fileformat=VCFv4.1
        '''
        s = ''
        for k in self.__spec:
            s += "##" +  "=".join([k, self.__spec[k]]) + "\n"
        return s
        
    def params_section(self):
        '''
        params_section()->str, returns VCF header of PARAM section
        e.g. ##PARAMS=<ID=is_qcfail,value=False>
        '''
        p = AlignmentConfig()
        params = ''
        for k in p.conf:
            params += "##PARAMS=" + ",".join(['<ID='+ k, 'Value='+str(p.conf[k])]) + '>\n'
        return params

    def info_section(self, id, number, value, desc):
        '''
        info_section(id, number, value, desc)->str, return INFO header
        e.g. ##INFO
        '''
        prefix = '##INFO='
        info_h = ''
        info_h += ','.join([prefix+ '<ID='+ id.upper(),
                            'Number=' + num,
                            'Type='+ value,
                            'Description='+ '"' + desc + '">\n'])
        return info_h
    
    def filter(self):
        return str('PASS')

    def header(self):
        return "{0}\t{1}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\n".format("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT")

    def vcf_stream(self):
        pass

        
if __name__ == '__main__':
    writer = VCFWriter("test")
    print writer.spec_section(),
    print writer.params_section(),
    print writer.header(),
    
