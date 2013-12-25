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
        return  len(reads)
        
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
    def mismatch_frequency(m=[], mis=[]):
        try:
            return len(mis) / (len(m) + len(mis))
        except ZeroDivisionError:
            return .0
            
    @staticmethod
    def a_to_g_frequency(a, g):
        try:
            return len(g) / (len(a) + len(g))
        except ZeroDivisionError:
            return .0
               
    @staticmethod
    def define_allele(base, ref=None):
        if base and ref:
            [_.upper() for _ in base]
            ref.upper()
        
        c = Counter(base)
        comm = c.most_common()

        __allele = {}
        for base in comm:
            if base[0] != ref:
                __allele.update({base[0]:base[1]})
        defined = ()
        for j in __allele:
            for k in __allele:
                # single alllele is found
                if j == k:
                    return tuple([j, __allele[j]])
                    
                # most common varinat with a allele type alone
                elif __allele[k] == __allele[j] and k != j:
                    return tuple(__allele.items())
                    
                # most common variant if has many allele
                elif __allele[k] != __allele[j] and k != j:
                    m = max(__allele[k], __allele[j])
                    if m == __allele[k]:
                        return tuple([k, __allele[k]])
                    elif m == __allele[j]:
                        return tuple([j, __allele[j]])
        else:
            # allele is not found
            return '.'

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
            if ref == 'A':
                ref_r = (ar)
                ref_f = (af)
                alt_f = (gf + tf + cf)
                alt_r = (gr + tr + cr)
            elif ref == 'T':
                ref_r = (tr)
                ref_f = (tf)
                alt_r = (gr + cr + ar)
                alt_f = (gf + cf + af)
            elif ref == 'G':
                ref_r = (gr)
                ref_f = (gf)
                alt_r = (cr + tr + ar)
                alt_f = (cf + cf + af)
            elif ref == 'C':
                ref_r = (cr)
                ref_f = (cf)
                alt_r = (ar + tr + gr)
                alt_f = (af + tf + gf)
            return tuple([ref_r, ref_f, alt_r, alt_f])
        else:
            raise ValueError(
                'Could not define the allele base {all_bases:s}, {chrom:s}, {pos:s}'
                .format(all_bases=all_bases, chrom=bam_chrom, pos=pos))
