from __future__ import division
from collections import Counter
import string
import re
import pysam

class AlignmentConfig(object):
    def __init__(self):
        self.conf = self.__set_default()
        
    def __set_default(self):
        __initialize = {
            'is_duplicate' : False,
            'is_unmapped' : False,
            'is_deletion' : False,
            'is_proper_pair' : True,
            'is_qcfail' : False,
            'is_secondary' : True,
            'mapq' : 25,
            'mate_is_reverse' : True,
            'mate_is_unmapped' : False,
            'base_qual' : 25,
        }
        return __initialize
        
    def set_filter(self, param, value):
        self.conf.update({param:value})
        
    def get_filter_value(self, param):
        return self.conf[param]
        
    def has_filter(self, filt):
        for k in self.conf:
            if filt == self.conf[k]:
                return True
            else:
                return False
        
    def print_all_params(self):
        for k in self.conf:
            print "[%s]:%s" % (k, self.conf[k])
            
        
class AlignmentStream(object):
    def __init__(self, samfile, fafile, chrom=None, start=None, end=None, one_based=None):
        bm = pysam.Samfile(samfile, 'rb')
        ft = pysam.Fastafile(fafile)
        
        self.samfile = bm
        self.fafile = ft
        self.chrom = chrom
        self.one_based = one_based
        (self.start, self.end) = self.__resolve_coords(start, end)

    def pileup_stream(self):
        mod_chr = ''
        if self.chrom.startswith('chr'):
            mod_chr = re.sub(r'^chr', '', self.chrom, 1)
        else:
            mod_chr = self.chrom
            
        for col in self.samfile.pileup(reference=mod_chr,
                                           start=self.start,
                                           end=self.end):
            if self.one_based:
                pos = col.pos + 1
            else:
                pos = col.pos
                
            prop_read = []
            for r in col.pileups:
                if r.alignment.is_proper_pair \
                    and not r.alignment.is_duplicate \
                    and not r.alignment.is_unmapped \
                    and not r.is_del:
                   prop_read.append(r)

            bam_chrom = self.samfile.getrname(col.tid)
            ref = ''
            if not bam_chrom.startswith('chr'):
                ref = self.fafile.fetch(reference=self.chrom, start=col.pos, end=col.pos+1)
            else:
                mod_chr = 'chr' + self.chrom
                ref = self.fafile.fetch(reference=mod_chr, start=col.pos, end=col.pos+1)
                
            if not ref:
                raise ValueError('No sequence content within [chr:%s, start:%s, end:%s]' % \
                                 (self.chrom, self.start, self.end))

            print mod_chr, self.chrom
                
            mismatches = [read for read in prop_read
                          if read.alignment.seq[read.qpos] != ref]
            if len(mismatches) > 1:
                
                # Mismatch base
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
                
                # returns per a base
                yield {'CHROM': chrom,
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


