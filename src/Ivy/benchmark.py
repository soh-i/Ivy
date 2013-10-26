from __future__ import division
import vcf
import os.path
import re
import ConfigParser
from collections import Counter
try:
    import utils
except:
    raise ImportError

    
class DarnedReader(object):
    '''
    DarnedReader generates the subset of DARNED db
    >>> darned_db = DarnedReader(sp='human_hg19', source='Brain', db='Path_to_Darned_DB')
    Returns list of subset of darned db
    '''
    def __init__(self, sp='', source=None, db_path=None):
        self.__sp = sp
        if source is not None:
            self.source = [k for k in source]
        else:
            self.source = None
        
        if db_path is None:
            conf_path = utils.find_app_root()
            conf = ConfigParser.RawConfigParser()
            conf.read(conf_path + '/data.ini')

            if conf.has_section('DARNED'):
                sp = conf.get('DARNED', self.__sp)
                self.__darned = {self.__sp : conf_path + sp}
                self.db = self.__generate_darned_set()
            else:
                raise (ConfigParser.NoSectionError, 'Invalid species name [%s] is given' % (self.__sp))

        elif db_path is not None:
            self.__darnd = {self.__sp:db_path}
            self.db = self.__generate_darned_set()

    def __str__(self):
        return 'Darned file [%s]' % (self.__darned[self.__sp])
        
    def __generate_darned_set(self):
        self.cnt = 0
        # Store selected records
        if self.__sp and self.source is not None:
            darned_list = []
            with open(self.__darned[self.__sp], 'r') as f:
                for line in f:
                    if not line.startswith('chrom'):
                        data = line.split(",")
                        chrom = data[0]
                        pos = data[1]
                        darned_source = data[8]
                        
                        if darned_source == self.source:
                            darned_list.append(chrom + ':' + pos + self.source)
                            self.cnt += 1
                return darned_list
                
        # Store all Darned records (default)
        elif self.__sp and self.source is None or self.source is 'all' or self.source is 'All':
            selected = []
            with open(self.__darned[self.__sp], 'r') as f:
                for line in f:
                    if not line.startswith('chrom'):
                        data = line.split(',')
                        chrom = data[0]
                        pos = data[1]
                        selected.append(chrom + ':' + pos + 'All')
                        self.cnt += 1
                return selected
                
        elif not self.__sp:
            raise (RuntimeError, 'Given species name[%s] is not valid' % (self.__sp))

    def sp(self):
        '''
        >>> darned_db.sp()
        -> str, given species name
        '''
        return self.__sp

    def path(self):
        '''
        >>> darned_db.path()
        -> str, absolute path to Darned database file
        '''
        return os.path.abspath(self.__darned[self.__sp])
        
    def db_name(self):
        '''
        >>> darned_db.db_name()
        -> str, Darned db name
        '''
        return os.path.basename(self.path())
                 
    def count(self):
        '''
        >>> darned_db.count(self)
        -> int, number of the Darned entories
        '''
        return self.cnt

        
class VCFReader(object):
    def __init__(self, filename):
        self.vcf = filename
        self.db = self.__generate_vcf_set()

    def __generate_vcf_set(self):
        '''
        generate_vcf_set(self) -> list, returns the accumulated vcf
        '''
        vcf_recs = []
        vcf_reader = vcf.Reader(open(self.vcf, 'r'))
        self.count = 0
        self.substitutions = Counter()
        
        for rec in vcf_reader:
            types = str(rec.REF) + '-to-' + 'or'.join([str(i) for i in rec.ALT])
            self.substitutions[types] += 1
            mod_chr = re.sub(r'^chr', '', rec.CHROM, 1)
            vcf_recs.append(mod_chr + ':' + str(rec.POS))
            self.count += 1
        return vcf_recs
        
    def cnt(self):
        '''
        count(self) -> int, entoties of the parsed vcf records
        '''
        return self.count

    def vcf_name(self):
        return os.path.basename(self.vcf)

    def editing_types(self):
        return self.substitutions
        
    def ag_count(self):
        return self.substitutions.get('A-to-G')
        
    def other_mutations_count(self):
        i = 0
        for k in self.substitutions:
            if not k == 'A-to-G':
                i += self.substitutions[k]
        return i
        
class Benchmark(object):
    def __init__(self, answer=None, predict=None):
        self.answer = set(answer)
        self.predict = set(predict)
        self.intersect = self.answer.intersection(self.predict)

    def precision(self):
        precision = 0
        try:
            precision = len(self.intersect)/len(self.predict)
        except ZeroDivisionError:
            pass
        finally:
            return precision
            
    def recall(self):
        recall = 0
        try:
            recall = len(self.intersect)/len(self.answer)
        except ZeroDivisionError:
            pass
        finally:
            return recall
            
    def f_measure(self):
        precision = self.precision()
        recall = self.recall()
        
        f = 0
        try:
            f = (2*recall*precision)/(recall+precision)
        except ZeroDivisionError:
            pass
        finally:
            return f

    def ag_enrichment(self):
        raise NotImplementedError
