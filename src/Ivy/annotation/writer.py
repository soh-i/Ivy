import sys
from Ivy.commandline.parse_ivy_opts import CommandLineParser
from Ivy.utils import AttrDict
from Ivy.version import __version__
import datetime
import os.path
import string
    

def to_xml(cls, value):
    return '##{cls:s}=<ID={_id:},Number={num:},Type={_type:},Description="{desc:}">\n'.format(
        cls=cls.upper(), _id=value.get("id"), num=value.get("num"), _type=value.get("type"), desc=value.get("desc"))
    
def header_column_name():
    return "{0}\t{1}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\n".format(
        "#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT")
    
    
class VCFInfoHeader(object):
    '''
    Args:
     vcf_info_data(dict): dic = {'id': hoge, 'num': 1, 'type': Integer, 'desc': piyo}
    
    Returns:
     #INFO=<ID=NS,Number=1,Type=Integer,Description="Number of Samples With Data>"
    '''
    
    def __init__(self):
        pass
    
    def build(self):
        s = self._NS()
        s += self._DP()
        s += self._DP4()
        s += self._AF()
        s += self._EF()
        s += self._MAPQ()
        s += self._BACQ()
        s += self._SB()
        s += self._PB()
        s += self._MA()
        s += self._MIS()
        return s

    def _NS(self):
        value = {'id': 'NS', 'num': 1, 'type': 'Integer', 'desc': 'Number of Samples With Data'}
        return to_xml("INFO", value)
        
    def _DP(self):
        value = {'id': 'DP', 'num': 1, 'type': 'Integer', 'desc': 'Total Depth'}
        return to_xml("INFO", value)

    def _DP4(self):
        value = {'id': 'DP4', 'num': 4, 'type': 'Integer',
                 'desc': 'ref-forward bases, ref-reverse, alt-forward and alt-reverse bases'}
        return to_xml("INFO", value)
        
    def _AF(self):
        value = {'id': 'AF', 'num': 1, 'type': 'Float', 'desc': 'Allele Frequency'}
        return to_xml("INFO", value)

    def _EF(self):
        value = {'id': 'EF', 'num': 1, 'type': 'Float', 'desc': 'Editing Frequency'}
        return to_xml("INFO", value)

    def _MAPQ(self):
        value = {'id': 'MAPQ', 'num': 1, 'type': 'Integer', 'desc': 'Average Mapping Quality'}
        return to_xml("INFO", value)

    def _BACQ(self):
        value = {'id': 'BACQ', 'num': 1, 'type': 'Integer', 'desc': 'Average Phread-scaled Base Call Quality'}
        return to_xml("INFO", value)
        
    def _SB(self):
        value = {'id': 'SB', 'num': 1, 'type': 'Float', 'desc': 'Strand Bias of P-value'}
        return to_xml("INFO", value)
        
    def _PB(self):
        value = {'id': 'PB', 'num': 1, 'type': 'Float', 'desc': 'Positional Bias of P-value'}
        return to_xml("INFO", value)
        
    def _MIS(self):
        value = {'id': 'MIS', 'num': 1, 'type': 'Integer', 'desc': 'Total Mismatch Reads'}
        return to_xml("INFO", value)
        
    def _MA(self):
        value = {'id': 'MA', 'num': 1, 'type': 'Integer', 'desc': 'Total Match Reads'}
        return to_xml("INFO", value)

class VCFMetaHeader(object):
    def __init__(self, __params):
        self.__spec = __params

    def build(self):
        header = self._add_date_spec()
        header += self._add_file_spec()
        header += self._add_params_spec()
        return header
        
    def _add_date_spec(self):
        d = datetime.datetime.today()
        _date_info = {
            'fileformat': 'VCFv4.1',
            'source': 'Ivy_' + __version__,
            'filedate': '{0}/{1}/{2} {3:0>2d}:{4:0>2d}:{5:0>2d}'.format(
                d.year, d.month, d.day, d.hour, d.minute, d.second)
        }
        s = str()
        for _ in _date_info:
            s += "##"+  "=".join(['{meta:s}'.format(
                meta=_), _date_info[_]]) + "\n"
        return s
        
    def _add_file_spec(self):
        files =  ["file://" + os.path.abspath(self.__spec.fasta),
                  "file://" +  os.path.abspath(self.__spec.r_bams)]
        if self.__spec.d_bams:
            files.append("file://" +  os.path.abspath(self.__spec.d_bams))
            
        return "".join(['##' + _ + '\n' for _ in files])
        
    def _add_params_spec(self):
        params = str()
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
    
        
class VCFHeaderBuilder(object):
    def __init__(self, params):
        self.params = params
        
    def build(self):
        meta = VCFMetaHeader(self.params)
        info = VCFInfoHeader()
        header = meta.build()
        header += info.build()
        header += header_column_name()
        return header

        
class VCFWriterDataLine(VCFInfoHeader):
    def __init__(self):
        pass

    def write_data_line(self, data):
        return('{chrom:}\t{pos:}\t{id_:}\t{ref:}\t{alt:}\t{qual:}\t{filt:}\t'.format(
            chrom=data.get('chrom'),
            pos=data.get('pos'),
            id_=data.get('id'),
            ref=data.get('ref'),
            alt=data.get('alt'),
            qual=data.get('qual'),
            filt=data.get('filt')))
        
    def _to_info(self, id_, val):
        return('{id_:}={val:};'.format(id_=id_, val=val))

    def write_info_line(self, data):
        s =  self._DP(data)
        s += self._DP4(data)
        s += self._AF(data)
        s += self._EF(data)
        s += self._MAPQ(data)
        s += self._BACQ(data)
        s += self._SB(data)
        s += self._PB(data)
        s += self._MA(data)
        s += self._MIS(data)
        return s
        
    def _DP(self, data):
        return self._to_info('DP', data.get('coverage'))

    def _DP4(self, data):
        return self._to_info('DP4', ",".join([str(_) for _ in data.get('dp4')]))

    def _AF(self, data):
        return self._to_info('AF', data.get('allele_freq'))

    def _EF(self, data):
        return self._to_info('EF', data.get('ag_freq'))

    def _MAPQ(self, data):
        return self._to_info('MAPQ', data.get('average_mapq'))

    def _BACQ(self, data):
        return self._to_info('BACQ', data.get('average_baq'))

    def _SB(self, data):
        return self._to_info('SB', data.get('strand_bias'))

    def _PB(self, data):
        return self._to_info('PB', data.get('positional_bias'))

    def _MIS(self, data):
        return self._to_info('MIS',  data.get('mismatches'))

    def _MA(self, data):
        return self._to_info('MA', data.get('matches'))
    
        


        
