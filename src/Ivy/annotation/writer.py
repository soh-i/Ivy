import sys
from alignment import AlignmentConfig

class VCFWriterException(Exception):
    pass

class VCFWriteHeader(object):
    def __init__(self):
        self.__spec = {
            'fileformat': 'VCFv4.1',
            'source': 'Ivy_v0.0.1',
            'reference': 'test.fasta',
            'species': 'Homo Sapiens',
            'samples': 'sample.bam',
        }

    def make_vcf_header(self):
        sys.stdout.write(self.__spec_section())
        sys.stdout.write(self.__params_section())
        sys.stdout.write(self.__header_name())

    def merge_filtering_param(self):
        '''
        merge_filtering_param() -> dict, returns merged filtering and spec dict
        '''
        p = AlignmentConfig()
        self.__spec.update({k: p.conf[k] for k in p.conf})
        return self.__spec

    def __spec_section(self):
        '''
        spec_section() -> str, return VCF header of basic informations
        e.g. ##fileformat=VCFv4.1
        '''
        s = ''
        for k in self.__spec:
            s += "##"+  "=".join([k, self.__spec[k]]) + "\n"
        return s
        
    def __params_section(self):
        '''
        params_section()->str, returns VCF header of PARAM section
        e.g. ##PARAMS=<ID=is_qcfail,value=False>
        '''
        p = AlignmentConfig()
        params = ''
        for k in p.conf:
            params += "##PARAMS="+ ",".join(['<ID='+ k, 'Value='+ str(p.conf[k])])+ '>\n'
        return params

    def __info_section(self, id, number, value, desc):
        '''
        info_section(id, number, value, desc)->str, return INFO header
        e.g. ##INFO
        '''
        prefix = '##INFO='
        info_h = ''
        info_h += ','.join([prefix+ '<ID='+ id.upper(),
                            'Number='+ num,
                            'Type='+ value,
                            'Description='+ '"'+ desc+ '">\n'])
        return info_h
    
    def __header_name(self):
        return "{0}\t{1}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\n".format("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT")

