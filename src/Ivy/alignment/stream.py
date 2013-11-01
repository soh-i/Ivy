from __future__ import division
from collections import Counter, namedtuple
import string
import re
import pysam
from Ivy.utils import ImutableDict

class AlignmentConfig(object):
    def __init__(self):
        self.params = self.__set_default()
        
    def __set_default(self):
       __params = {
           'is_duplicate': False,
           'is_unmapped': False,
           'is_deletion': False,
           'is_proper_pair': True,
           'is_qcfail': False,
           'is_secondary': True,
           'mapq': 25,
           'mate_is_reverse': True,
           'mate_is_unmapped': False,
           'base_qual': 25,
           'edit_ratio': 0.1,
           'edit_base_c': 10,
           'mutation_type_c': 1
       }
       self.__imutable_conf = ImutableDict(__params)
       return self.__imutable_conf

       
class AlignmentStream(object):
    def __init__(self, samfile, fafile, chrom=None, start=None, end=None, one_based=True):
        __bm = pysam.Samfile(samfile, 'rb')
        __ft = pysam.Fastafile(fafile)
        
        self.samfile = __bm
        self.fafile = __ft
        self.one_based = one_based
        (self.start, self.end) = self.__resolve_coords(start, end, one_based)

        if not chrom.startswith('chr'):
            self.chrom = 'chr' + chrom
        else: self.chrom = chrom

        debug = False
        if debug:
            # info. for loaded samfile
            print self.samfile.filename
            print self.samfile.lengths
            print self.samfile.mapped
            print self.samfile.nreferences
            print self.samfile.references
            print self.samfile.unmapped
            # info. for fasta
            print self.fafile.filename

    def alignment_prepare(self):
        raise NotImplementedError
        
    def __sort(self):
        pass

    def __index(self):
        pass

    def __faidx(self):
        pass
            
    def pileup_stream(self):
        for col in self.samfile.pileup(reference=self.chrom,
                                       start=self.start,
                                       end=self.end):
            bam_chrom = self.samfile.getrname(col.tid)
            if self.one_based:
                pos = col.pos + 1
            else:
                pos = col.pos
            #pos = col.pos + 1 if self.one_based else col.pos
            
            ref = self.fafile.fetch(reference=bam_chrom, start=col.pos, end=col.pos+1).upper()
            reads = col.pileups
            
            # Raw reads (no filterings through)
            raw_reads = [_ for _ in reads]
            raw_mismatches = [_ for _ in raw_reads if _.alignment.seq[_.qpos] != ref]
            raw_matches = [_ for _ in raw_reads if _.alignment.seq[_.qpos] == ref]

            # Has proper_pair and without deletion
            prop_nodel_reads = [_ for _ in reads if not _.is_del and _.alignment.is_proper_pair]
            prop_nodel_mismatchs = [_ for _ in prop_nodel_reads if _.alignment.seq[_.qpos] != ref]
            prop_nodel_matches = [_ for _ in prop_nodel_reads if _.alignment.seq[_.qpos] == ref]
            
            # Has NO deletion
            nodel_reads = [_ for _ in reads if not _.is_del]
            nodel_mismatches = [_ for _ in nodel_reads if _.alignment.seq[_.qpos] != ref]
            nodel_matches = [_ for _ in nodel_reads if _.alignment.seq[_.qpos] == ref]

            # Has proper_pair alone
            prop_reads = [_ for _ in reads if _.alignment.is_proper_pair]
            prop_mismatches =  [_ for _ in prop_reads if _.alignment.seq[_.qpos] != ref]
            prop_matches =  [_ for _ in prop_reads if _.alignment.seq[_.qpos] == ref]
            
            # Has deletions alone
            del_reads = [_ for _ in reads if not _.is_del]
            del_prop_reads = [_ for _ in reads if not _.is_del]

            # Has insertion alone
            ins_reads = [_ for _ in reads if _.is_del > 0]
            ins_prop_reads = [_ for _ in reads if _.is_del > 0]
            
            filt_reads = []
            for _ in col.pileups:
                if _.alignment.is_proper_pair \
                   and not _.alignment.is_secondary \
                   and not _.alignment.is_qcfail \
                   and not _.alignment.is_duplicate \
                   and not _.alignment.is_unmapped \
                   and not _.is_del:
                    filt_reads.append(_)
                    
            filt_mismatches = [_ for _ in filt_reads if _.alignment.seq[_.qpos] != ref]
            filt_matches = [_ for _ in filt_reads if _.alignment.seq[_.qpos] == ref]

            
            if not ref:
                raise ValueError('No seq. content within [chr:%s, start:%s, end:%s]' % \
                                 (self.chrom, self.start, self.end))
                        
            A = [_ for _ in prop_reads if _.alignment.seq[_.qpos] == 'A']
            C = [_ for _ in prop_reads if _.alignment.seq[_.qpos] == 'C']
            T = [_ for _ in prop_reads if _.alignment.seq[_.qpos] == 'T']
            G = [_ for _ in prop_reads if _.alignment.seq[_.qpos] == 'G']
            N = [_ for _ in prop_reads if _.alignment.seq[_.qpos] == 'N']
            
            G_r = [_.alignment.is_reverse for _ in G
                   if _.alignment.is_reverse]
            g_f = [_.alignment.is_reverse for _ in G
                   if not _.alignment.is_reverse]
            A_r = [_.alignment.is_reverse for _ in A
                   if _.alignment.is_reverse]
            a_f = [_.alignment.is_reverse for _ in A
                   if not _.alignment.is_reverse]
            T_r = [_.alignment.is_reverse for _ in T
                   if _.alignment.is_reverse]
            t_f = [_.alignment.is_reverse for _ in T
                   if not _.alignment.is_reverse]
            C_r = [_.alignment.is_reverse for _ in C
                   if _.alignment.is_reverse]
            c_f = [_.alignment.is_reverse for _ in C
                   if not _.alignment.is_reverse]
            N_r = [_.alignment.is_reverse for _ in N
                   if _.alignment.is_reverse]
            n_f = [_.alignment.is_reverse for _ in N
                   if not _.alignment.is_reverse]

            mutation_type = {'A': len(A), 'T': len(T), 'G': len(G), 'C': len(C)}
            
                        
            debug = False
            if debug:
                coverage = A_r+ a_f+ T_r+ t_f+ G_r+ g_f+ C_r+ c_f + N_r + n_f
                #print [_.alignment.seq[_.qpos] for _ in G]
                print '[A:%s,%s] [T:%s,%s] [G:%s,%s] [C:%s,%s]' % (A_r, a_f, T_r, t_f, G_r, g_f, C_r, c_f)
                print 'Coverage:%d' % (coverage)
                
            try:
                pass
                #allele_freq = mismatch_c / len(raw_reads)
                #ag_freq = (mc_G + mc_g) / (mc_G + mc_g + mc_A + mc_a)
            except ZeroDivisionError:
                pass
                
            yield {
                'CHROM': self.chrom,
                'POS': pos,
                'REF': ref,
                'raw_coverage': len(raw_reads),
                'prop_coverage': len(prop_reads),
                'prop_nodel_coverage': len(prop_nodel_reads),
                'nodel_coverage': len(nodel_reads),
                'mismatches': len(raw_mismatches),
                'matches': len(raw_matches),
                'cov': len(raw_mismatches) + len(raw_matches),
                'types': mutation_type,
                'A': len(A),
                'T': len(T),
                'C': len(C),
                'G': len(G),
                'N': len(N),
            }
    
    def __resolve_coords(self, start, end, is_one_based):
        if is_one_based:
            if start is not None:
                start = start -1
            else:
                None
            if end is not None:
                end = end -1
            else:
                None
        return start, end

    def average_baq(self, string):
        return [ord(s)-33 for s in string]

    def define_allele(self, base, ref=None):
        c = Counter(base)
        comm = c.most_common()
        allele = []
        for base in comm:
            if base[0] != ref:
                allele.append(base)
        return allele


def define_allele(base, ref=None):
    if base and ref:
        [_.upper() for _ in base]
        ref.upper()
        
    c = Counter(base)
    comm = c.most_common()
    #return [base for base in comm if base != ref]
    
    allele = {}
    for base in comm:
        if base[0] != ref:
            allele.update({base[0]:base[1]})
            
    for j in allele:
        for k in allele:
            # Has a allele alone
            if allele[k] == allele[j] and k != j:
                return tuple(allele)
                
            # Has many allele
            elif allele[k] != allele[j] and k != j:
                return k,j
        
if __name__ == '__main__':
    a = ['A', 'T', 'C', 'G']
    b = ['C', 'G', 'G', 'G', 'A', 'A', 'A']
    c = ['A', 'T', 'C', 'G']
    print a, 'r:A',
    print define_allele(a, ref='A') #=> C, T, Gx

    print b, 'r:G',
    print define_allele(b, ref='G') #=> A

    print c, 'r:A',
    print define_allele(c, ref='A')
        
        
class AlignmentStreamMerger(object):
    def __init__(self, rna, dna):
        self.rna = rna
        self.dna = dna

    def merge_streaming(self):
        dna_stream = AlignmentStream(conf)
        rna_stream = AlignmentStream(conf)


