from __future__ import division
import vcf

class FileReader(object):
    def __init__(self, file):
        self.vcf = file

    def genrate_vcf_set(self):
        '''
        returns set of given VCF
        '''
        pass

    def generate_darned_set(self):
        '''
        returns set of Darned db
        '''
        pass

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
        return (2*recall*precision)/(recall+precision)
        
