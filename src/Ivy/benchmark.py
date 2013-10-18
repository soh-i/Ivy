from __future__ import division
import vcf
import os.path
import ConfigParser

class DarnedReader(object):
    def __init__(self, sp=''):
        '''
        args: sp='' 
        '''
        self.__sp = sp
        self.__darned = {
            'human':'../data/DARNED_hg19.csv',
            'fly':'',
            'mice':''
        }

    def generate_darned_set(self):
        '''
        args: species
        return: set of Darned list, chr:pos
        '''
        darned_list = []
        self.count = 0
        if self.__sp == 'Human' or self.__sp == 'human':
            with open(self.__darned[self.__sp], 'r') as f:
                for line in f:
                    if not line.startswith('chrom'):
                        data = line.split(",")
                        chr = data[0]
                        pos = data[1]
                        darned_list.append(chr+':'+pos)
                        self.count += 1
            return darned_list
        else:
            raise (RuntimeError, 'Given species name[%s] is not valid' % (self.__sp))

    def sp(self):
        return self.__sp

    def path(self):
        return os.path.abspath(self.__darned[self.__sp])

    def db_name(self):
        return os.path.basename(self.path())
                 
    def count(self):
        return self.count
        

class VCFReader(object):
    def __init__(self, file):
        self.vcf = file

    def generate_vcf_set(self):
        '''
        args: None
        returns: set of VCF record list
        '''
        vcf_recs = []
        self.count = 0
        vcf_reader = vcf.Reader(open(self.vcf, 'r'))
        for rec in vcf_reader:
            vcf_recs.append(rec.CHROM + ':' + str(rec.POS))
            self.count += 0
        return vcf_recs

    def count(self):
        return self.count

    def vcf_name(self):
        return os.path.basename(self.vcf)

        
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
