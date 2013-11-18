import sys
from Ivy.parse_opt import CommandLineParser
from Ivy.utils import AttrDict
import datetime

class VCFWriteHeader(object):
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
        AttrDict.show(self.__spec)
        
    def make_vcf_header(self):
        sys.stdout.write(self.__spec_section())
        sys.stdout.write(self.__params_section())
        sys.stdout.write(self.__header_name())

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
            params += '##PARAMS=' + ','.join([
                '<ID={id},Value={val},Filter_class={filt}>\n'.format(
                    id=_, val=self.__spec.basic_filter._data[_], filt='basic_filter')])

        # ext filter group
        for _ in self.__spec.ext_filter._data:
            params += '##PARAMS=' + ','.join([
                '<ID={id},Value={val},Filter_class={filt}>\n'.format(
                    id=_, val=self.__spec.ext_filter._data[_], filt='ext_filter')])

        # stat filter group
        for _ in self.__spec.stat_filter._data:
            params += '##PARAMS=' + ','.join([
                '<ID={id},Value={val},Filter_class={filt}>\n'.format(
                    id=_, val=self.__spec.stat_filter._data[_], filt='stat_filter')])
            
        # sample filter group
        for _ in self.__spec.sample_filter._data:
            params += '##PARAMS=' + ','.join([
                '<ID={id},Value={val},Filter_class={filt}>\n'.format(
                    id=_, val=self.__spec.sample_filter._data[_], filt='sample_filter')])
             
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
        return "{0}\t{1}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\n".format(
            "#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT")

