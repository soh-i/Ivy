from __future__ import division
from collections import Counter, namedtuple
import os.path
import string
import re
import math
import pprint
import logging
import warnings

__program__ = 'stream'
__author__ = 'Soh Ishiguro <yukke@g-language.org>'
__license__ = ''
__status__ = 'development'


class AlignmentReadsStats(object):
    '''
    Utility class provides methods to streaming/filtering reads processing as staticmethods,
    This class can NOT to handle the Pysam object.
    
    Example:
     >>> import AlignmentUtils
     >>> AlignmentUtils.mismatch_frequency(match, mismatch)
    '''
    
    @staticmethod
    def quals_in_pos(reads):
        return [ord(_.alignment.qual[_.qpos])-33 for _ in reads]
        
    @staticmethod
    def reads_coverage(reads):
        return len(reads)
        
    @staticmethod
    def average_base_quality(reads):
        q_pos = AlignmentReadsStats.quals_in_pos(reads)
        try:
            return math.ceil(sum(q_pos)/len(q_pos))
        except ZeroDivisionError:
            return .0

    @staticmethod
    def average_mapq(reads):
        mapqs_in_pos = [_.alignment.mapq for _ in reads]
        try:
            return math.ceil(sum(mapqs_in_pos) / len(mapqs_in_pos))
        except ZeroDivisionError:
            return .0

    @staticmethod
    def mismatch_frequency(m=None, mis=None):
        if isinstance(m, list) and isinstance(mis, list):
            try:
                return len(mis) / (len(m) + len(mis))
            except ZeroDivisionError:
                return .0
            
    @staticmethod
    def a_to_g_frequency(a=None, g=None):
        if isinstance(a, list) and isinstance(g, list):
            try:
                return len(g) / (len(a + g))
            except ZeroDivisionError:
                return .0
               
    @staticmethod
    #def define_allele(base, ref=None):
    #    if base and ref:
    #        [_.upper() for _ in base]
    #        ref.upper()
    #    
    #    c = Counter(base)
    #    comm = c.most_common()
    # 
    #    __allele = {}
    #    for base in comm:
    #        if base[0] != ref:
    #            __allele.update({base[0]:base[1]})
    # 
    #    for j in __allele:
    #        for k in __allele:
    #            # single alllele is found
    #            if j == k:
    #                return tuple([j, __allele[j]])
    #                
    #            # most common varinat with a allele type alone
    #            elif __allele[k] == __allele[j] and k != j:
    #                return tuple(__allele.items())
    #                
    #            # most common variant if has many allele
    #            elif __allele[k] != __allele[j] and k != j:
    #                m = max(__allele[k], __allele[j])
    #                if m == __allele[k]:
    #                    return tuple([k, __allele[k]])
    #                elif m == __allele[j]:
    #                    return tuple([j, __allele[j]])
    #    else:
    #        # allele is not found
    #        return '.'
    
    def define_allele(base, ref=None):
        found_allele = {}
        data = Counter(base)
        for i in data.most_common():
            b_type = i[0]
            count = i[1]
            if b_type != ref:
                found_allele.update({b_type: count})
        if len(found_allele):
            # Mismatch is found
            max_val= max(found_allele.values())
        else:
            # Mismatch base is not found
            return '.'
            
        result = []
        for value in found_allele.values():
            if value == max_val:
                for itm in found_allele.items():
                    if itm[1] == value:
                        result.append(itm)
                break
        return result
        
    @staticmethod
    def compute_dp4(ref=None, af=None, ar=None, tf=None, tr=None, gf=None, gr=None, cf=None, cr=None):
        '''
        Args:
         number of base with specific direction(int)
        Returns:
         dp4(tuple) tuple has elements in
          Number of 1) forward ref alleles
                    2) reverse ref
                    3) forward non-ref
                    4) reverse non-ref alleles
        '''
        if af + ar + gf + gr + cf + cr + tf + tr != 0:
            ref_r = 0
            ref_f = 0
            alt_r = 0
            alt_f = 0
            dp4 = 0
            if ref == 'A':
                # REF: A
                # ALT: T, G, C
                ref_f = (af)
                ref_r = (ar)
                alt_f = (gf + tf + cf)
                alt_r = (gr + tr + cr)
                dp4 = tuple([ref_f, ref_r, alt_f, alt_r])
                
            elif ref == 'T':
                # REF: T
                # ALT: A, G, C
                ref_f = (tf)
                ref_r = (tr)
                alt_f = (gf + cf + af)
                alt_r = (gr + cr + ar)
                dp4 = tuple([ref_f, ref_r, alt_f, alt_r])
                
            elif ref == 'G':
                # REF: G
                # ALT: A, T, C
                ref_f = (gf)
                ref_r = (gr)
                alt_f = (cf + tf + af)
                alt_r = (cr + tr + ar)
                dp4 = tuple([ref_f, ref_r, alt_f, alt_r])
                
            elif ref == 'C':
                # REF: C
                # ALT: A, T, G
                ref_f = (cf)
                ref_r = (cr)
                alt_f = (af + tf + gf)
                alt_r = (ar + tr + gr)
                dp4 = tuple([ref_f, ref_r, alt_f, alt_r])
            return dp4
            
        else:
            raise ValueError(
                'Could not define the allele base {all_bases:s}, {chrom:s}, {pos:s}'
                .format(all_bases=all_bases, chrom=bam_chrom, pos=pos))

if __name__ == '__main__':
    dp4 = AlignmentReadsStats.compute_dp4(ref="A",
                                          af=11, ar=1,
                                          tf=1, tr=0,
                                          gf=2, gr=9,
                                          cf=9, cr=1)
    ma = [_ for _ in range(1, 598)]
    mis = [_ for _ in range(1,31)]
    freq = AlignmentReadsStats.mismatch_frequency(mis=mis, m=ma)
    print dp4
    print freq
    
    base = ['A', 'T', 'T', 'T', 'T', 'G', 'G', 'G', 'A', 'A']
    allele = AlignmentReadsStats.define_allele(base, ref='T')
    print base
    print "REF: T"
    print allele
    
    
    
