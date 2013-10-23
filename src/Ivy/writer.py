import collections
from alignment import AlignmentConfig

class VCFWriter(object):
    def __init__(self):
        
        self.spec = {'fileformat': 'VCFv4.1',
                     'source': 'Ivy_v0.0.1',
                     'reference': 'test.fasta'
                   }
    
    def add_to_spec(self):
        p = AlignmentConfig()
        self.spec.update({k: p.conf[k] for k in p.conf})
        return self.spec

    def spec_section(self):
        '''
        spec_section() -> str, return VCF header of basic informations
        e.g. ##fileformat=VCFv4.1
        '''
        s = ''
        for k in self.spec:
            s += "##" +  "=".join([k, self.spec[k]]) + "\n"
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
        info_h += ','.join([prefix+ '<ID='+ id,
                            'Number=' + num,
                            'Type='+ value,
                            'Description='+ '"' + desc + '">\n'])
        return info_h
    
    def format(self):
        for i in range(10):
            yield i

    def filter(self):
        return str('PASS')

        
if __name__ == '__main__':
    writer = VCFWriter()
    print writer.spec_section(),

    
    
    

