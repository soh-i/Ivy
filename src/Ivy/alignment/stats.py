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
    def compute_dp4(ref, A_f, A_r, T_f, T_r, G_f, G_r, C_f, C_r):
        #if len(alt): TODO:  here is bug # TypeError: object of type 'NoneType' has no len()
        if True: # TODO: set any condition(s)
            ref_r = 0
            ref_f = 0
            alt_r = 0
            alt_f = 0
            
            if ref == 'A':
                ref_r = (A_r)
                ref_f = (A_f)
                alt_r = (G_f+C_f+T_f)
                alt_f = (G_r+C_r+T_r)
            elif ref == 'T':
                ref_r = (T_r)
                ref_f = (T_f)
                alt_r = (G_r+C_r+A_r)
                alt_f = (G_f+C_f+A_f)
            elif ref == 'G':
                ref_r = (G_r)
                ref_f = (G_f)
                alt_r = (C_r+T_r+A_r)
                alt_f = (C_f+C_f+C_r)
            elif ref == 'C':
                ref_r = (C_r)
                ref_f = (C_f)
                alt_r = (A_r+T_r+G_r)
                alt_f = (A_f+T_f+G_f)
            return tuple([ref_r, ref_f, alt_r, alt_f])

        else:
            raise RuntimeError(
                'Could not able to define the allele base {all_bases:s}, {chrom:s}, {pos:s}'
                .format(all_bases=all_bases, chrom=bam_chrom, pos=pos))


