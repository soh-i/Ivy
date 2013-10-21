from __future__ import division
import vcf
import os.path
import ConfigParser
from collections import Counter
try:
    import utils
except:
    raise ImportError

class DarnedReader(object):
    '''
    DardReader class
    '''
    def __init__(self, sp=''):
        '''
        __init__(self, sp=) -> instance object of the DarnedReader
        '''
        self.__sp = sp
        conf_path = utils.find_app_root()
        conf = ConfigParser.RawConfigParser()
        conf.read(conf_path + '/data.ini')

        if conf.has_section('DARNED'):
            sp = conf.get('DARNED', self.__sp)
            print sp
            self.__darned = {self.__sp : conf_path + sp}
            self.db = self.__generate_darned_set()
        else:
            raise (ConfigParser.NoSectionError, 'Invalid species name [%s] is given' % (self.__sp))
            
    def __str__(self):
        return 'Darned file [%s]' % (self.__darned[self.__sp])
        
    def __generate_darned_set(self):
        '''
        generate_darned_set(self) -> list, returns the accumulated DARNED db
        '''
        darned_list = []
        self.cnt = 0
        if self.__sp:
            with open(self.__darned[self.__sp], 'r') as f:
                for line in f:
                    if not line.startswith('chrom'):
                        data = line.split(",")
                        chr = data[0]
                        pos = data[1]
                        darned_list.append(chr + ':' + pos)
                        self.cnt += 1
            return darned_list
        else:
            raise (RuntimeError, 'Given species name[%s] is not valid' % (self.__sp))

    def sp(self):
        '''
        sp(self) -> str, given species name
        '''
        return self.__sp

    def path(self):
        '''
        path(self) -> str, absolute path to Darned database file
        '''
        return os.path.abspath(self.__darned[self.__sp])

    def db_name(self):
        '''
        db_name(self) -> str, Darned db name
        '''
        return os.path.basename(self.path())
                 
    def count(self):
        '''
        count(self) -> int, number of the Darned entories
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
            vcf_recs.append(rec.CHROM + ':' + str(rec.POS))
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
        pass


