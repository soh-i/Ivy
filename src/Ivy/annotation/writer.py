import sys
from Ivy.commandline.parse_ivy_opts import CommandLineParser
from Ivy.utils import AttrDict
import datetime
import os.path
from functools import wraps

def to_xml(value):
    if not isinstance(value, dict):
        return False
        
    def _xml(function):
        @wraps(function)
        
        def __xml(*args, **kw):
            result = function(*args, **kw)
            return '##INFO=<ID={id:},Number={num:},Type={type:},Description="{desc:}">'.format(
                id=value.get("id"), num=value.get("num"), type=value.get("type"), desc=value.get("desc"))
        return __xml
    return _xml

    
class VCFWriter(object):
    def __init__(self):
        pass
        
class VCFWriterInfoHeader(VCFWriter):
    '''
    Args:
     vcf_info_data(dict): dic = {'id': hoge, 'num': 1, 'type': Integer, 'desc': piyo}
    
    Returns:
     #INFO=<ID=NS,Number=1,Type=Integer,Description="Number of Samples With Data>"
    '''
    
    def __init__(self):
        VCFWriter.__init__(self)

    def _to_xml(self, value):
        return '##INFO=<ID={id:},Number={num:},Type={type:},Description="{desc:}">'.format(
            id=value.get("id"), num=value.get("num"), type=value.get("type"), desc=value.get("desc"))

    def write_info_header(self):
        print self._NS()
        print self._DP()
        print self._DP4()
        print self._AF()
        print self._EF()
        print self._MAPQ()
        print self._BACQ()
        print self._SB()
        print self._PB()
        print self._MA()
        print self._MIS()

    def _NS(self):
        value = {'id': 'NS', 'num': 1, 'type': 'Integer', 'desc': 'Number of Samples With Data'}
        return self._to_xml(value)
        
    def _DP(self):
        value = {'id': 'DP', 'num': 1, 'type': 'Integer', 'desc': 'Total Depth'}
        return self._to_xml(value)

    def _DP4(self):
        value = {'id': 'DP4', 'num': 4, 'type': 'Integer', 'desc': 'ref-forward bases, ref-reverse, alt-forward and alt-reverse bases'}
        return self._to_xml(value)
        
    def _AF(self):
        value = {'id': 'AF', 'num': 1, 'type': 'Float', 'desc': 'Allele Frequency'}
        return self._to_xml(value)

    def _EF(self):
        value = {'id': 'EF', 'num': 1, 'type': 'Float', 'desc': 'Editing Frequency'}
        return self._to_xml(value)

    def _MAPQ(self):
        value = {'id': 'MAPQ', 'num': 1, 'type': 'Integer', 'desc': 'Average Mapping Quality'}
        return self._to_xml(value)

    def _BACQ(self):
        value = {'id': 'BACQ', 'num': 1, 'type': 'Integer', 'desc': 'Average Phread-scaled Base Call Quality'}
        return self._to_xml(value)
        
    def _SB(self):
        value = {'id': 'SB', 'num': 1, 'type': 'Float', 'desc': 'Strand Bias of P-value'}
        return self._to_xml(value)
        
    def _PB(self):
        value = {'id': 'PB', 'num': 1, 'type': 'Float', 'desc': 'Positional Bias of P-value'}
        return self._to_xml(value)
        
    def _MIS(self):
        value = {'id': 'MIS', 'num': 1, 'type': 'Integer', 'desc': 'Total Mismatch Reads'}
        return self._to_xml(value)
        
    def _MA(self):
        value = {'id': 'MA', 'num': 1, 'type': 'Integer', 'desc': 'Total Match Reads'}
        return self._to_xml(value)

class VCFWriterInfoDataLine(VCFWriterInfoHeader):
    pass

        
class VCFWriterDataLine(VCFWriterInfoHeader):
    def __init__(self, data):
        self.data = data

    def chrom(self):
        pass

    def pos(self):
        pass

    def id(self):
        pass

    def ref(self):
        pass

    def alt(self):
        pass

    def qual(self):
        pass

    def filter(self):
        pass

    def info(self):
        pass
    
    def format(self):
        pass
        
class VCFWriterMetaInformationLine(VCFWriter):
    def __init__(self):
        pass

        
class VCFWriteHeader(VCFWriter):
    def __init__(self, __params):
        d = datetime.datetime.today()
        _info = {
            'fileformat': 'VCFv4.1',
            'source': 'Ivy_v0.0.1',
            'filedate': '{0}/{1}/{2} {3}:{4}:{5}'.format(
                d.year, d.month, d.day, d.hour, d.minute, d.second)
        }
        __params._vcf_meta = _info
        self.__spec = __params
        self.__spec._vcf_meta.filename = "file://" + os.path.abspath(self.__spec.fasta)
        self.__spec._vcf_meta.bam = "file://" +  os.path.abspath(self.__spec.r_bams)
        
    def make_vcf_header(self):
        sys.stdout.write(self.__spec_section())
        sys.stdout.write(self.__params_section())
        #sys.stdout.write(self.header_column_name())

    def merge_filtering_param(self):
        self.__spec.update({k: p.conf[k] for k in p.conf})
        return self.__spec

    def __spec_section(self):
        s = ''
        for _ in self.__spec._vcf_meta._data:
            s += "##"+  "=".join(['{meta:s}'.format(
                meta=_), self.__spec._vcf_meta[_]]) + "\n"
        return s
            
    def __params_section(self):
        params = ''
        # basic filter group
        for _ in self.__spec.basic_filter._data:
            params += '##IVY_PARAMS=' + ','.join([
                '<ID={id},Value={val},Class={filt}>\n'.format(
                    id=_, val=self.__spec.basic_filter._data[_], filt='basic_filter')])

        # ext filter group
        for _ in self.__spec.ext_filter._data:
            params += '##IVY_PARAMS=' + ','.join([
                '<ID={id},Value={val},Class={filt}>\n'.format(
                    id=_, val=self.__spec.ext_filter._data[_], filt='ext_filter')])

        # stat filter group
        for _ in self.__spec.stat_filter._data:
            params += '##IVY_PARAMS=' + ','.join([
                '<ID={id},Value={val},Class={filt}>\n'.format(
                    id=_, val=self.__spec.stat_filter._data[_], filt='stat_filter')])
            
        # sample filter group
        for _ in self.__spec.sample_filter._data:
            params += '##IVY_PARAMS=' + ','.join([
                '<ID={id},Value={val},Class={filt}>\n'.format(
                    id=_, val=self.__spec.sample_filter._data[_], filt='sample_filter')])
        return params

    def __info_section(self, id, number, value, desc):
        prefix = '##INFO='
        info_h = ''
        info_h += ','.join([prefix+ '<ID='+ id.upper(),
                            'Number='+ num,
                            'Type='+ value,
                            'Description='+ '"'+ desc+ '">\n'])
        return info_h
    
    def header_column_name(self):
        sys.stdout.write("{0}\t{1}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\n".format(
            "#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT"))

