from __future__ import division
from collections import Counter, namedtuple
import os.path
import string
import re
import pysam
from Ivy.utils import ImutableDict

__program__ = 'stream'
__author__ = 'Soh Ishiguro <yukke@g-language.org>'
__license__ = ''
__status__ = 'development'


class AlignmentConfig(object):
    def __init__(self, params_from_cl):
        #self.params = self.__set_default()
        
        if isinstance(params_from_cl, dict):
            self.cl_params = ImutableDict(params_from_cl)
        else:
            raise ValueError, ('Dict is only acceptable as command-line argument')
        
    #def __set_default(self):
    #   __params = {
    #       'is_duplicate': False,
    #       'is_unmapped': False,
    #       'is_deletion': False,
    #       'is_proper_pair': True,
    #       'is_qcfail': False,
    #       'is_secondary': True,
    #       'mapq': 25,
    #       'mate_is_reverse': True,
    #       'mate_is_unmapped': False,
    #       'base_qual': 25,
    #       'edit_ratio': 0.1,
    #       'edit_base_c': 10,
    #       'mutation_type_c': 1
    #   }
    #   self.__imutable_conf = ImutableDict(__params)
    #   return self.__imutable_conf
       
    def logger(self):
        # TODO: logging for used params with value
        pass

    def config_varidator(self):
        # TODO: to varidate each params/values
        pass

        
class AlignmentPreparation(object):
    def __init__(self):
        # TODO: AlignmentStream.__init__ move into here.
        pass
        
    def alignment_prepare(self):
        raise NotImplementedError
    
    def __sort(self):
        if not os.path.isfile(bamfile):
            try:
                pysam.sort(self.samfile, self.samfile + 'sorted')
                sort_log = pysam.sort.getMessage()
                return True
            except:
                raise RuntimeError()
        else:
            print "already sorted"
            return False

    def __index(self):
        if not os.path.isfile(samfile + '.index.bam'):
            try:
                pysam.index(self.samfile)
                return True
            except:
                raise RuntimeError()
        else:
            print "already indexed"
            return False

    def __faidx(self):
        if not os.path.isfile(fafile + '.fai'):
            try:
                pysam.faidx(self.fafile)
                return True
            except:
                raise RuntimeError()
        else:
            print "already exist"
            return False 

    def __merge_bams(self, bams=[]):
        for _ in bams:
            if not os.path.isfile(_):
                raise RuntimeError()
        try:
            pysam.merge([_ for _ in bams])
            return True
        except:
            raise RuntimeError()

def __resolve_chrom_name(bam, fa):
    bam_obj = pysam.Samfile(bam, 'rb')
    fa_obj = pysam.Fastafile(fa)
    

class AlignmentStream(object):
    def __init__(self, config):
        __bm = pysam.Samfile(config['r_bams'], 'rb', check_header=True, check_sq=True)
        __ft = pysam.Fastafile(config['fasta'])
        
        self.samfile = __bm
        self.fafile = __ft
        self.one_based = config['one_based']
        (self.start, self.end) = self.__resolve_coords(config['region']['start'], config['region']['end'], self.one_based)

        if not config['region']['chrom'].startswith('chr'):
            self.chrom = 'chr' + config['region']['chrom']
        else: self.chrom = config['region']['chrom']
        
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
        
    def pileup_stream(self):
        for col in self.samfile.pileup(reference=self.chrom,
                                       start=self.start,
                                       end=self.end,
                                       ):
            
            bam_chrom = self.samfile.getrname(col.tid)
            if self.one_based:
                pos = col.pos + 1
            else:
                pos = col.pos
            
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

            # Has insertion alone
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

            if not ref:
                # TODO: resolve difference name in fasta and bam
                raise ValueError('No seq. content within [chr:%s, start:%s, end:%s]' % \
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
        return start, end

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
                        
if __name__ == '__main__':
    conf = AlignmentConfig()
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
    ## 
    ## 
    #   
       
class AlignmentStreamMerger(object):
    def __init__(self, rna, dna):
        self.rna = rna
        self.dna = dna

    def merge_streaming(self):
        dna_stream = AlignmentStream(conf)
        rna_stream = AlignmentStream(conf)


