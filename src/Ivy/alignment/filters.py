import sys
import fisher

__program__ = 'filter'
__author__ = 'Soh Ishiguro <yukke@g-language.org>'
__license__ = ''
__status__ = 'development'


def strand_bias_filter(mr, mf, mir, mif):
    p = fisher.pvalue(mr, mf, mir, mif)
    return p
    
        
        

















    
