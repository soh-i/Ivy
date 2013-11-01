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
       
if __name__ == '__main__':
    align = AlignmentConfig()
    
    align.params.replace("mapq", 33)
    print align.params["mapq"]
    print align.params
        
class AlignmentStream(object):
    def __init__(self, samfile, fafile, chrom=None, start=None, end=None, one_based=True):
        __bm = pysam.Samfile(samfile, 'rb')
        __ft = pysam.Fastafile(fafile)
        
        self.samfile = __bm
        self.fafile = __ft
        self.one_based = one_based
        (self.start, self.end) = self.__resolve_coords(start, end)

        if not chrom.startswith('chr'):
            self.chrom = 'chr' + chrom
        else: self.chrom = chrom

        debug = True
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
            
    def pileup_stream(self):
        for col in self.samfile.pileup(reference=self.chrom,
                                       start=self.start,
                                       end=self.end):
            if self.one_based:
                pos = col.pos + 1
            else: pos = col.pos
                
            prop_read = []
            for r in col.pileups:
                #if r.alignment.is_proper_pair \
                #    and not r.alignment.is_duplicate \
                #    and not r.alignment.is_unmapped \
                #    and not r.is_del:
                prop_read.append(r)
            
            bam_chrom = self.samfile.getrname(col.tid)
            ref = self.fafile.fetch(reference=bam_chrom, start=col.pos, end=col.pos+1)
            
            if not ref:
                raise ValueError('No seq. content within [chr:%s, start:%s, end:%s]' % \
                                 (self.chrom, self.start, self.end))
                        
            mismatches = [_ for _ in prop_read
                          if _.alignment.seq[_.qpos] != ref]
            if len(mismatches) > 1:
                
                # Mismatch basen
                A = [read for read in prop_read
                     if read.alignment.seq[read.qpos] == 'A']
                a = [read for read in prop_read
                     if read.alignment.seq[read.qpos] == 'a']
                C = [read for read in prop_read
                     if read.alignment.seq[read.qpos] == 'C']
                c = [read for read in prop_read
                     if read.alignment.seq[read.qpos] == 'c']
                T = [read for read in prop_read
                     if read.alignment.seq[read.qpos] == 'T']
                t = [read for read in prop_read
                     if read.alignment.seq[read.qpos] == 't']
                G = [read for read in prop_read
                     if read.alignment.seq[read.qpos] == 'G']
                g = [read for read in prop_read
                     if read.alignment.seq[read.qpos] == 'g']
                N = [read for read in prop_read
                     if read.alignment.seq[read.qpos] == 'N' \
                     or read.alignment.seq[read.qpos] == 'n']

                # separeted by read directions
                debug = True
                if debug:
                    G_r = len([_.alignment.is_reverse for _ in G
                               if _.alignment.is_reverse])
                    g_f = len([_.alignment.is_reverse for _ in G
                               if not _.alignment.is_reverse])
                    
                    A_r = len([_.alignment.is_reverse for _ in A
                               if _.alignment.is_reverse])
                    a_f = len([_.alignment.is_reverse for _ in A
                               if not _.alignment.is_reverse])
                    
                    T_r = len([_.alignment.is_reverse for _ in T
                               if _.alignment.is_reverse])
                    t_f = len([_.alignment.is_reverse for _ in T
                               if not _.alignment.is_reverse])
                    
                    C_r = len([_.alignment.is_reverse for _ in C
                               if _.alignment.is_reverse])
                    c_f = len([_.alignment.is_reverse for _ in C
                               if not _.alignment.is_reverse])
                    
                    N_r = len([_.alignment.is_reverse for _ in N
                               if _.alignment.is_reverse])
                    n_f = len([_.alignment.is_reverse for _ in N
                               if not _.alignment.is_reverse])
                    coverage = A_r+ a_f+ T_r+ t_f+ G_r+ g_f+ C_r+ c_f + N_r + n_f
                    
                    print [_.alignment.seq[_.qpos] for _ in g] # is not working
                    print [_.alignment.seq[_.qpos] for _ in G]
                    
                    print '[A:%s,%s] [T:%s,%s] [G:%s,%s] [C:%s,%s]' % (A_r, a_f, T_r, t_f, G_r, g_f, C_r, c_f)
                    print 'Coverage:%d' % (coverage)
                    
                mc_A = len(A)
                mc_a = len(a)
                mc_T = len(T)
                mc_t = len(t)
                mc_G = len(G)
                mc_g = len(g)
                mc_C = len(C)
                mc_c = len(c)
                mc_N = len(N)
                
                Ac = mc_A + mc_a
                Tc = mc_T + mc_t
                Gc = mc_G + mc_g
                Cc = mc_C + mc_c
                base = {'A': Ac, 'T': Tc, 'G': Gc, 'C': Cc}
                alleles = self.define_allele(base, ref=ref)
                
                forward_allel_c = (mc_A + mc_T + mc_C + mc_G)
                reverse_allel_c = (mc_a + mc_t + mc_c + mc_g)
                depth = len(prop_read)

                mismatch_c = sum([base[1] for base in alleles])
                mismatch_freq = '{0:.2f}'.format(mismatch_c/depth)
                
                allele_freq = 0
                ag_freq = 0
                try:
                    allele_freq = mismatch_c / depth
                    ag_freq = (mc_G + mc_g) / (mc_G + mc_g + mc_A + mc_a)
                except ZeroDivisionError:
                    pass
                
                mapq = r.alignment.mapq
                ave_baq = '{0:.2f}'.format(sum(self.average_baq(r.alignment.seq))/depth)
                #print self.average_baq(r.alignment.seq)
                
                # returns per a base
                yield {
                    'CHROM': self.chrom,
                    'POS': pos,
                    'REF': ref,
                    'ALT': ",".join([(b[0]) for b in alleles]),
                    'ID': 'ID',
                    'FORMAT': '.',
                    'INFO': '.',
                    'FILTER': '.',
                    'coverage': len(prop_read),
                    'mismatches': mismatch_c,
                    'mismatch_freq': mismatch_freq,
                    'matches': depth,
                    'mapq': mapq,
                    'QUAL': ave_baq,
                    'forward_allel_count': forward_allel_c,
                    'reverse_allel_count': reverse_allel_c,
                    'Af':mc_A,
                    'Ar':mc_a,
                    'Cf':mc_C,
                    'Cr':mc_c,
                    'Tf':mc_T,
                    'Tr':mc_t,
                    'Gf':mc_G,
                    'Gr':mc_g,
                    'N': mc_N,
                }
                
    def __resolve_coords(self, start, end):
        if self.one_based:
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

        
class AlignmentStreamMerger(object):
    def __init__(self, rna, dna):
        self.rna = rna
        self.dna = dna

    def merge_streaming(self):
        dna_stream = AlignmentStream(conf)
        rna_stream = AlignmentStream(conf)


