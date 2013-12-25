import sys
import fisher
from scipy import stats, mean

__program__ = 'filter'
__author__ = 'Soh Ishiguro <yukke@g-language.org>'
__license__ = ''
__status__ = 'development'

def strand_bias_filter(m=None, mis=None):
    mis_r = [_.alignment.seq[_.qpos] for _ in mis if _.alignment.is_reverse]
    mis_f = [_.alignment.seq[_.qpos] for _ in mis if not _.alignment.is_reverse]
    ma_r =  [_.alignment.seq[_.qpos] for _ in m if _.alignment.is_reverse]
    ma_f =  [_.alignment.seq[_.qpos] for _ in m if not _.alignment.is_reverse]
    p = fisher.pvalue(len(ma_r), len(ma_f), len(mis_r), len(mis_f))
    return p.left_tail

def positional_bias_filter(m=None, mis=None):
    mismatch_pos = [_.alignment.qlen for _ in mis]
    match_pos = [_.alignment.qlen for _ in m]
    try:
        p = stats.ttest_ind(match_pos, mismatch_pos)[1]
    except ZeroDivisionError:
        p = .0
    return p
        
    
















    
