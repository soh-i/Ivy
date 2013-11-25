import sys
import fisher

__program__ = 'filter'
__author__ = 'Soh Ishiguro <yukke@g-language.org>'
__license__ = ''
__status__ = 'development'

def strand_bias_filter(match, mismatch):
    mis_r = [_.alignment.seq[_.qpos] for _ in mismatch if _.alignment.is_reverse]
    mis_f = [_.alignment.seq[_.qpos] for _ in mismatch if not _.alignment.is_reverse]
    ma_r =  [_.alignment.seq[_.qpos] for _ in match if _.alignment.is_reverse]
    ma_f =  [_.alignment.seq[_.qpos] for _ in match if not _.alignment.is_reverse]
    p = fisher.pvalue(len(ma_r), len(ma_f), len(mis_r), len(mis_f))
    return p.left_tail
        
def positional_bias_filter(match, mismatch):
    pass
    
















    
