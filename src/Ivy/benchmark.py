from __future__ import division
import vcf

class DarnedReader(object):
    def __init__(self, sp=''):
        '''
        args: sp='' 
        '''
        self.__sp = sp
        self.__darned = {
            'human':'/home/soh.i/benchmark/data/human/DARNED_hg19.csv',
            'fly':'',
            'mice':''
        }

    def generate_darned_set(self):
        '''
        args: species
        return: set of Darned list, chr:pos
        '''
        darned_list = []
        if self.__sp == 'Human' or self.__sp == 'human':
            with open(self.__darned[self.__sp], 'r') as f:
                for line in f:
                    if not line.startswith('chrom'):
                        data = line.split(",")
                        chr = data[0]
                        pos = data[1]
                        darned_list.append(chr+':'+pos)
            return darned_list
        else:
            raise (RuntimeError, 'Given species name[%s] is not valid' % (self.__sp))

            
class VCFReader(object):
    def __init__(self, file):
        self.vcf = file

    def generate_vcf_set(self):
        '''
        args: None
        returns: set of VCF record list
        '''
        vcf_recs = []
        vcf_reader = vcf.Reader(open(self.vcf, 'r'))
        for rec in vcf_reader:
            vcf_recs.append(rec.CHROM + ':' + str(rec.POS))
        return vcf_recs

        
class Benchmark(object):
    def __init__(self, answer=None, predict=None):
        self.answer = set(answer)
        self.predict = set(predict)
        self.intersect = self.answer.intersection(self.predict)

    def precision(self):
        if not len(self.predict) == 0 or not len(self.answer) == 0:
            precision = len(self.intersect)/len(self.predict)
            return precision
        else:
            raise ValueError('Error: set object has 0 element')
            
    def recall(self):
        if not len(self.predict) == 0 or not len(self.answer) == 0:
            recall = len(self.intersect)/len(self.answer)
            return recall
        else:
            raise ValueError('Error: set object has 0 element')
            
    def f_measure(self):
        precision = self.precision()
        recall    = self.recall()
        f = 0
        try:
            f = (2*recall*precision)/(recall+precision)
        except ZeroDivisionError:
            pass
            
        finally:
            return f
