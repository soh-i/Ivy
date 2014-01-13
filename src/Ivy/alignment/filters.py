import sys

__program__ = 'filter'
__author__ = 'Soh Ishiguro <yukke@g-language.org>'
__license__ = ''
__status__ = 'development'

def strand_bias_filter(m=None, mis=None):
    import fisher
    
    mis_r = [_.alignment.seq[_.qpos] for _ in mis if _.alignment.is_reverse]
    mis_f = [_.alignment.seq[_.qpos] for _ in mis if not _.alignment.is_reverse]
    ma_r =  [_.alignment.seq[_.qpos] for _ in m if _.alignment.is_reverse]
    ma_f =  [_.alignment.seq[_.qpos] for _ in m if not _.alignment.is_reverse]
    p = fisher.pvalue(len(ma_r), len(ma_f), len(mis_r), len(mis_f))
    return p.left_tail

def positional_bias_filter(m=None, mis=None):
    from scipy import stats, mean
    
    mismatch_pos = [_.alignment.qlen for _ in mis]
    match_pos = [_.alignment.qlen for _ in m]
    try:
        p = stats.ttest_ind(match_pos, mismatch_pos)[1]
    except ZeroDivisionError:
        p = .0
    return p

def base_call_bias(m=None, mis=None):
    pass


def gtest(f_obs, f_exp=None, ddof=0):
    import numpy as np
    from scipy.stats import chisqprob, chisquare
    f_obs = np.asarray(f_obs, 'f')
    k = f_obs.shape[0]
    f_exp = np.array([np.sum(f_obs, axis=0) / float(k)] * k, 'f') \
            if f_exp is None \
               else np.asarray(f_exp, 'f')
    g = 2 * np.add.reduce(f_obs * np.log(f_obs / f_exp))
    return g, chisqprob(g, k - 1 - ddof)

if __name__ == '__main__':
    print  gtest([10,2], [2,10])

    
    












