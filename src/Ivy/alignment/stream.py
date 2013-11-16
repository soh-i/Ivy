from __future__ import division
from collections import Counter, namedtuple
import os.path
import string
import re
import pysam
from Ivy.utils import die

__program__ = 'stream'
__author__ = 'Soh Ishiguro <yukke@g-language.org>'
__license__ = ''
__status__ = 'development'


class AlignmentStream(object):
    def __init__(self, __params):
        if hasattr(__params, 'AttrDict'):
            self.params = __params
        else:
            raise TypeError("Given param %s is %s class, not 'AttrDic' class" %
                            (__params, __params.__class__.__name__))

        __bm = pysam.Samfile(self.params.r_bams, 'rb', check_header=True, check_sq=True)
        __ft = pysam.Fastafile(self.params.fasta)
        
        self.samfile = __bm
        self.fafile = __ft
        self.one_based = self.params.one_based

        ### Resolve to explore specified region or not
        # explore all region, set None
        #die(self.params)
        if self.params.region.all_flag == 'All':
            print "ALL"
            #self.params.region.start = 1
            #self.params.region.end =1
            #self.params.region.chrom = 1

        # explore specified region
        elif self.params.region.chrom and self.params.region.start and self.params.region.end:
            (self.params.region.start, self.params.region.end) = self.__resolve_coords(
                self.params.region.start,
                self.params.region.end,
                self.params.one_based)
            if not self.params.region.chrom.startswith('chr'):
                self.params.region.chrom = 'chr' + self.params.region.chrom
            else: self.params.region.chrom = self.params.region.chrom
        else:
            #raise ValueError("chrom or pos of start/end is not set")
            pass

        debug = False
        if debug:
            # info. for loaded samfile
            print "### info. for samfile object from given Bam header @SQ ###"
            print "Sam file: %s" % self.samfile.filename
            print "lengths: %s" % [_ for _ in self.samfile.lengths]
            print "mapped %d: " % self.samfile.mapped
            print "N_references: %s" % self.samfile.nreferences
            print "references: %s" % [_ for _ in self.samfile.references]
            print "unmapped: %s" % self.samfile.unmapped
            
            # info. for fasta
            print "### info. for fasfile object ###"
            print "filename: %s" % self.fafile.filename

        if _is_same_chromosome_name(bam=self.params.r_bams, fa=self.params.fasta):
            pass
        else:
            raise RuntimeError("invalid chrom name")
            
    def pileup_stream(self):
        print self.params.region.chrom
        print self.params.region.start
        print self.params.region.end
        
        for col in self.samfile.pileup(reference=self.params.region.chrom,
                                       start=self.params.region.start,
                                       end=self.params.region.end
                                       ):
            
            bam_chrom = self.samfile.getrname(col.tid)
            if self.params.one_based:
                pos = col.pos + 1
            else:
                pos = col.pos
                
            #print self.fafile.fetch(reference=21, start=int(col.pos), end=int(col.pos)+1)
            ref = self.fafile.fetch(reference=bam_chrom, start=col.pos,
                                    end=col.pos+1).upper()

            
            reads = col.pileups
            
            # Raw reads (no filterings through)
            #raw_reads = [_ for _ in reads]
            #raw_mismatches = [_ for _ in raw_reads if _.alignment.seq[_.qpos] != ref]
            #raw_matches = [_ for _ in raw_reads if _.alignment.seq[_.qpos] == ref]
            # 
            ## Has proper_pair and without deletion
            #prop_nodel_reads = [_ for _ in reads if not _.is_del and _.alignment.is_proper_pair]
            #prop_nodel_mismatchs = [_ for _ in prop_nodel_reads if _.alignment.seq[_.qpos] != ref]
            #prop_nodel_matches = [_ for _ in prop_nodel_reads if _.alignment.seq[_.qpos] == ref]
            # 
            ## Has NO deletion
            #nodel_reads = [_ for _ in reads if not _.is_del]
            #nodel_mismatches = [_ for _ in nodel_reads if _.alignment.seq[_.qpos] != ref]
            #nodel_matches = [_ for _ in nodel_reads if _.alignment.seq[_.qpos] == ref]
            # 
            ## Has proper_pair alone
            #prop_reads = [_ for _ in reads if _.alignment.is_proper_pair]
            #prop_mismatches =  [_ for _ in prop_reads if _.alignment.seq[_.qpos] != ref]
            #prop_matches =  [_ for _ in prop_reads if _.alignment.seq[_.qpos] == ref]
            # 
            ## Has deletions alone
            #del_reads = [_ for _ in reads if not _.is_del]
            #del_prop_reads = [_ for _ in reads if not _.is_del]

            # Has insertion alonep
            # TODO: fixt to print pysam object directory
            ins_reads = [_ for _ in reads if _.is_del > 0]
            ins_prop_reads = [_ for _ in reads if _.is_del > 0]
            del_reads = [_ for _ in reads if _.is_del < 0]
            del_prop_reads = [_ for _ in reads if _.is_del < 0]
            
            filt_reads = []
            for _ in col.pileups:
                if _.alignment.is_proper_pair \
                   and not _.alignment.is_secondary:
                   #and not _.alignment.is_qcfail \
                   #and not _.alignment.is_duplicate \
                   #and not _.alignment.is_unmapped \
                   #and not _.is_del:
                    filt_reads.append(_)
                    
            filt_mismatches = [_ for _ in filt_reads if _.alignment.seq[_.qpos] != ref]
            filt_matches = [_ for _ in filt_reads if _.alignment.seq[_.qpos] == ref]

            #print self.chrom, self.start, self.end
            if not ref:
                # TODO: resolve difference name in fasta and bam
                raise ValueError('No seq. content within [chr:%s, start:%s, end:%s], maybe different name of fasta and bam' % \
                                 (self.chrom, self.start, self.end))

            # array in read object per base types
            A = [_ for _ in filt_reads if _.alignment.seq[_.qpos] == 'A']
            C = [_ for _ in filt_reads if _.alignment.seq[_.qpos] == 'C']
            T = [_ for _ in filt_reads if _.alignment.seq[_.qpos] == 'T']
            G = [_ for _ in filt_reads if _.alignment.seq[_.qpos] == 'G']
            N = [_ for _ in filt_reads if _.alignment.seq[_.qpos] == 'N']

            # base string for 4 nucleotide types
            Gb =  [_.alignment.seq[_.qpos] for _ in G]
            Ab =  [_.alignment.seq[_.qpos] for _ in A]
            Tb =  [_.alignment.seq[_.qpos] for _ in T]
            Cb =  [_.alignment.seq[_.qpos] for _ in C]
            
            # strand informations
            G_r = [_.alignment.is_reverse for _ in G
                   if _.alignment.is_reverse].count(True)
            G_f = [_.alignment.is_reverse for _ in G
                   if not _.alignment.is_reverse].count(False)
            
            A_r = [_.alignment.is_reverse for _ in A
                   if _.alignment.is_reverse].count(True)
            A_f = [_.alignment.is_reverse for _ in A
                   
                   if not _.alignment.is_reverse].count(False)
            T_r = [_.alignment.is_reverse for _ in T
                   if _.alignment.is_reverse].count(True)
            T_f = [_.alignment.is_reverse for _ in T
                   if not _.alignment.is_reverse].count(False)
            
            C_r = [_.alignment.seq[_.qpos] for _ in C
                   if _.alignment.is_reverse].count(True)
            C_f = [_.alignment.seq[_.qpos] for _ in C
                   if not _.alignment.is_reverse].count(False)
            
            N_r = [_.alignment.is_reverse for _ in N
                   if _.alignment.is_reverse].count(True)
            N_f = [_.alignment.is_reverse for _ in N
                   if not _.alignment.is_reverse].count(False)

            mutation_type = ({'A': len(A), 'T': len(T), 'G': len(G),
                              'C': len(C), 'N': len(N)})

            Ac = [_.alignment.seq[_.qpos] for _ in A].count('A')
            Tc = [_.alignment.seq[_.qpos] for _ in T].count('T')
            Gc = [_.alignment.seq[_.qpos] for _ in G].count('G')
            Cc = [_.alignment.seq[_.qpos] for _ in C].count('C')
            Nc = [_.alignment.seq[_.qpos] for _ in C].count('N')
            coverage = Ac + Tc + Gc + Cc + Nc
            
            _all_base = Ab + Gb + Cb + Tb
            alt = self.define_allele(_all_base, ref=ref)
            
            # compute DP4 collumn
            # TODO: to write unittest is needed!
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
                elif ref == 'N':
                    ref_r = (N_r)
                    ref_f = (N_f)
                    alt_r = (A_r+T_r+G_r+C_r)
                    alt_f = (A_f+T_f+G_f+C_f)
                dp4 = tuple([ref_r, ref_f, alt_r, alt_f])
                
            else:
                raise RuntimeError(
                    'Could not able to define the allele base %s, chr[%s], pos[%s]'
                    % (all_bases, bam_chrom, pos))
            
            debug = False
            if debug:
                coverage = A_r+ a_f+ T_r+ t_f+ G_r+ g_f+ C_r+ c_f + N_r + n_f
                #print [_.alignment.seq[_.qpos] for _ in G]
                print '[A:%s,%s] [T:%s,%s] [G:%s,%s] [C:%s,%s]' \
                    % (A_r, a_f, T_r, t_f, G_r, g_f, C_r, c_f)
                print 'Coverage:%d' % (coverage)
            
            try:
                allele_ratio= len(filt_mismatches) / (len(filt_mismatches) + len(filt_matches))
                ag_ratio = len(G) / (len(G) + len(A))
            except ZeroDivisionError:
                allele_ratio = float(0)
                ag_ratio = float(0)

            yield {
                'chrom': bam_chrom,
                'pos': pos,
                'ref': ref,
                'alt': alt,
                'coverage': len(filt_reads),
                'mismatches': len(filt_mismatches),
                'matches': len(filt_matches),
                'cov': coverage,
                'mismach_ratio': allele_ratio,
                'ag_ratio': ag_ratio,
                'types': mutation_type,
                'Ac': len(A),
                'Tc': len(T),
                'Cc': len(C),
                'Gc': len(G),
                'Nc': len(N),
                'Gr': (G_r),
                'Gf': (G_f),
                'Cr': (C_r),
                'Cf': (C_f),
                'Tf': (T_f),
                'Tr': (T_r),
                'Af': (A_f),
                'Ar': (A_r),
                'Nr': (N_r),
                'Nf': (N_f),
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
        return int(start), int(end)

    def average_baq(self, string):
        return [ord(s)-33 for s in string]

    @classmethod
    def define_allele(self, base, ref=None):
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
            
   
class AlignmentStreamMerger(object):
    def __init__(self, rna, dna):
        self.rna = rna
        self.dna = dna

    def merge_streaming(self):
        dna_stream = AlignmentStream(conf)
        rna_stream = AlignmentStream(conf)


def _is_same_chromosome_name(bam=None, fa=None):
    __bam = pysam.Samfile(os.path.abspath(bam), 'rb')
    __fa = pysam.Fastafile(os.path.abspath(fa))
    bam_references = __bam.references
    fa_filename = __fa.filename
    fa_dx_filename = fa_filename + '.fai'
    
    if os.path.isfile(fa_dx_filename):
        for bam_chr in bam_references:
            if any([fai_chr == bam_chr for fai_chr in _parse_faidx(fa_dx_filename)]):
                return True
            else:
                return False
    else:
        raise RuntimeError("%s of faidx file is not found" % (fa_dx_filename))
        
def _parse_faidx(filename):
    fasta_chrom_name = []
    with open(filename, 'r') as fh:
        for row in fh:
            data = row.split('\t')
            fasta_chrom_name.append(data[0])
    return fasta_chrom_name

    
if __name__ == '__main__':
    #print _is_same_chromosome_name(bam="../data/testREDItools/dna.bam", fa="../data/testREDItools/reference.fa")
    test = AlignmentStream()
    print dir(test)

    #conf = AlignmentConfig()
    #a = ['A', 'T', 'C', 'G']
    #b = ['C', 'G', 'G', 'G', 'A', 'A', 'A']
    #c = ['A', 'T', 'C', 'G']
    ##d = ['A', 'A', 'A', 'T', 'T', 'T', 'T', 'C', 'C', 'G']
    # 
    #d = ['A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'T']
    #base = ['G', 'G', 'G', 'G', 'G', 'G', 'G', 'G', 'G', 'T']
    # 
    ##print a, 'r:A',
    #ref = ['A', 'T', 'G', 'C', 'N']
    #print base
    #for _ in ref:
    #    print "ref:", _,
    #    print AlignmentStream.define_allele(base, ref=_)
    # 
    ###print b, 'r:G',
    ##print define_allele(b, ref='G') #=> A
    ## 
    ###print c, 'r:A',
    ##print define_allele(c, ref='A')
    ## 
    ###print d, 'r:G',
    ##print define_allele(d, ref='G')
    
