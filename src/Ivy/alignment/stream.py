from __future__ import division
from collections import Counter, namedtuple
import os.path
import string
import re
import math
import pprint
import logging
import warnings
import pysam

from Ivy.utils import die, AttrDict, IvyLogger
from Ivy.alignment.filters import strand_bias_filter
from Ivy.alignment.stats import AlignmentReadsStats

__program__ = 'stream'
__author__ = 'Soh Ishiguro <yukke@g-language.org>'
__license__ = ''
__status__ = 'development'

DEBUG = False

            
class FilteredAlignmentReadsGenerator(object):
    '''
    AlignmentReadsFilter provides reads filter methods
    
    Args:
     Pysam.pileup object (common of all methods)
    Returns:
     reads(list), matches(list), mismatches(list)
    '''

    def reads_filter_by_all_params(self, pileup_obj, ref_base):
        '''
        Required attributes:
         proper pair reads
         NOT duplicated reads
         NOT unmapped reads
         NOT deletion reads
         NOT fail to quality check reads
        '''
                
        _reads = [_ for _ in pileup_obj
                  if (_.alignment.is_proper_pair
                      and not _.alignment.is_qcfail
                      and not _.alignment.is_duplicate
                      and not _.alignment.is_unmapped
                      and not _.is_del)]
        _matches = [_ for _ in _reads if _.alignment.seq[_.qpos] == ref_base]
        _mismatches = [_ for _ in _reads if _.alignment.seq[_.qpos] != ref_base]
        return _reads, _matches, _mismatches
            
    def reads_filter_without_pp(self, pileup_obj, ref_base):
        '''
        Required attributes:
         proper pair reads
         NOT fail to quality check reads
         NOT duplicated reads
         NOT unmapped reads
         NOT deletion reads
        '''
        
        _reads = [_ for _ in pileup_obj
                  if (_.alignment.is_proper_pair
                      and not _.alignment.is_qcfail
                      and not _.alignment.is_duplicate
                      and not _.alignment.is_unmapped
                      and not _.is_del)]
        _mismatches = [_ for _ in _reads if _.alignment.seq[_.qpos] != ref_base]
        _matches = [_ for _ in _reads if _.alignment.seq[_.qpos] == ref_base]
        return _reads, _matches, _mismatches
        
    def reads_allow_duplication(self, pileup_obj, ref_base):
        '''
        Required attributes:
         proper pair reads
         NOT failt to quality check reads
         NOT unmapped reads
         NOT deletion reads
        '''
        
        _reads = [_ for _ in col.pileups
                  if (_.alignment.is_proper_pair
                      and not _.alignment.is_qcfail
                      and not _.alignment.is_unmapped
                      and not _.is_del)]
        _matches = [_ for _ in _reads if _.alignment.seq[_.qpos] == ref_base]
        _mismatches = [_ for _ in _reads if _.alignment.seq[_.qpos] != ref_base]
        return _reads, _matches, _mismatches
            
    def reads_without_filter(self, pileup_obj, ref_base):
        '''All reads are passed through under this filter'''
        
        _reads = [_ for _ in pileup_obj if (not _.alignment.is_unmapped)]
        _mismatches = [_ for _ in _reads if _.alignment.seq[_.qpos] != ref_base]
        _matches = [_ for _ in _reads if _.alignment.seq[_.qpos] == ref_base]
        return _reads, _matches, _mismatches
        
    def reads_allow_insertion(self):
        print "Use --rm-insertion-reads is recommended"
        return None
            
    def reads_allow_deletion(self, pileup_obj, ref_base):
        print "Use --rm-deletion-reads is recommended"
        return None

        
class AlignmentStream(FilteredAlignmentReadsGenerator):
    def __init__(self, __params):
        '''
        Initialize for pileup bam files to explore RDD sites
        
        Args:
         Alignment parameters wrapped by AttrDic object
        
        Raises:
         TypeError:
         RuntimeError:
        
        Attributes:
         samfile:
         fasfile:
         one_based:
         params:
        '''
        
        ig = IvyLogger()
        self.logger = logging.getLogger(type(self).__name__)
        
        if hasattr(__params, 'AttrDict'):
            self.params = self.__add_preset(__params)
        else:
            raise TypeError("Given param {prm:s} is '{cls:s}' class, not 'AttrDic' class"
                            .format(prm=__params, cls=__params.__class__.__name__))

        self.samfile = pysam.Samfile(self.params.r_bams, 'rb', check_header=True, check_sq=True)
        self.fafile = pysam.Fastafile(self.params.fasta)
        self.one_based = self.params.one_based

        ### Resolve to explore specified region or not
        # explore all region if all_flag == 1
        if self.params.region.all_flag == 1:
            self.params.region.start = None
            self.params.region.end = None
            self.params.region.chrom = None

        # explore specified region if all_flag == 0
        elif self.params.region.all_flag == 0:
            if (self.params.region.chrom
                and self.params.region.start and self.params.region.end):
                (self.params.region.start, self.params.region.end) = \
                self.__resolve_coords(
                    self.params.region.start,
                    self.params.region.end,
                    self.params.one_based)
                
                if not self.params.region.chrom.startswith('chr'):
                    self.params.region.chrom = 'chr' + self.params.region.chrom
                else: self.params.region.chrom = self.params.region.chrom
            else:
                # explore all region if self.params.region.* is None
                pass
        if  _is_same_chromosome_name(bam=self.params.r_bams, fa=self.params.fasta):
            pass
        else:
            raise ValueError("Invalid chromosome name")

        if self.params.verbose:
            self.logger.debug(AttrDict.show(self.params))
            
        if DEBUG:
            _fasta_info()
            _sam_info()

    def __configuration(self):
        self.params.debug('Configuration...')
        pass
    
    def __add_preset(self, __p):
        '''
        Preset params add to the AttrDic attribute
        Examples:
         >>> __add_preset(AttrDic object)
        Returns:
         AttrDic object
        '''
        
        __p.preset_params._skip_N = True
        __p.preset_params._skip_has_many_allele = True
        return __p
            
    def pileup_stream(self):
        if self.params.verbose:
            self.logger.debug("Start pileup bam file...")
        
        for col in self.samfile.pileup(reference=self.params.region.chrom,
                                       start=self.params.region.start,
                                       end=self.params.region.end
                                       ):
            bam_chrom = self.samfile.getrname(col.tid)
            if self.params.one_based:
                pos = col.pos + 1
            else:
                pos = col.pos
            ref_base= self.fafile.fetch(reference=bam_chrom,
                                        start=col.pos,
                                        end=col.pos+1).upper()
            if not ref_base:
                # TODO: resolve difference name in fasta and bam
                raise ValueError(
                    'No sequence content within {chrom:s}, {start:s}, {end:s}'.format(
                        chrom=self.chrom, start=self.start, end=self.end))
            elif ref_base == 'N' or ref_base == 'n':
                continue
            
            #####################################
            ### Loading alignment with params ###
            #####################################
            params_debug = False
            #filter reads with all params
            if (self.params.basic_filter.rm_duplicated
                and self.params.basic_filter.rm_deletion
                and self.params.basic_filter.rm_insertion):
                if params_debug:
                    self.params.show(self.params.basic_filter)
                    raise SystemExit("Called methods: {0:s}".format(self.reads_filter_by_all_params.__name__))
                    
                passed_reads, passed_matches, passed_mismatches = self.reads_filter_by_all_params(col.pileup, ref_base)
                
            # allow duplicated containing reads
            elif (not self.params.basic_filter.rm_duplicated
                  and self.params.basic_filter.rm_deletion
                  and self.params.basic_filter.rm_insertion):
                if params_debug:
                    self.params.show(self.params.basic_filter)
                    raise SystemExit("Called methods: {0:s}".format(self.reads_allow_duplication.__name__))
                passed_reads, passed_matches, passed_mismaches = self.reads_allow_duplication(col.pileup, ref_base)
                
            # allow deletions containing reads
            elif (not self.params.basic_filter.rm_deletion
                  and self.params.basic_filter.rm_insertion
                  and self.params.basic_filter.rm_duplicated):
                if params_debug:
                    self.params.show(self.params.basic_filter)
                    raise SystemExit("Called methods: {0:s}".format(self.reads_allow_deletion.__name__))
                passed_reads, passed_matches, passed_mismatches = self.reads_allow_deletion(col.pileup, ref_base)
                
            # allow insertion containing reads
            elif (not self.params.basic_filter.rm_insertion
                  and self.params.basic_filter.rm_deletion
                  and self.params.basic_filterf.rm_duplicated):
                if params_debug:
                    self.params.show(self.params.basic_filter)
                    raise SystemExit("Called methods: {0:s}".format(self.params.basic_filter, self.reads_allow_insertion.__name__))
                passed_reads, passed_mathces, passed_mismatches = self.reads_allow_insertion(col.pileup, ref_base)

            # no filter
            else:
                if params_debug:
                    self.params.show(self.params.basic_filter)
                    raise SystemExit("Called methods: {0:s}".format(self.reads_without_filter.__name__))
                passed_reads, passed_matches, passed_mismatches = self.reads_without_filter(col.pileups, ref_base)

            ##############################
            ### Basic filters in reads ###
            ##############################
            # --min-rna-baq
            quals_in_pos = AlignmentReadsStats.quals_in_pos(passed_reads)
            average_baq = AlignmentReadsStats.average_base_quality(passed_reads)
            
            # --min-rna-cov
            coverage = AlignmentReadsStats.reads_coverage(passed_reads)

            # --min-rna-mapq
            average_mapq = AlignmentReadsStats.average_mapq(passed_reads)

            A, T, G, C = self.retrieve_reads_each_base_type(passed_reads)
            
            if (self.params.basic_filter.min_rna_cov <= coverage
                and self.params.basic_filter.min_rna_mapq <= average_mapq
                and self.params.basic_filter.min_baq_rna <= average_baq):
                
                # --min-mis-frequency
                allele_freq= AlignmentReadsStats.mismatch_frequency(passed_matches, passed_mismatches)
                ag_freq = AlignmentReadsStats.a_to_g_frequency(A, G)
                
                # --num-allow-type
                mutation_type = {}
                if len(A) > 0 and ref_base != 'A':
                    mutation_type.update({'A': len(A)})
                elif len(T) > 0 and ref_base != 'T':
                    mutation_type.update({'T': len(T)})
                elif len(G) > 0 and ref_base != 'G':
                    mutation_type.update({'G': len(G)})
                elif len(C) > 0 and ref_base != 'C':
                    mutation_type.update({'C': len(C)})
                    
                if (len(mutation_type) <= self.params.basic_filter.num_type
                    and len(mutation_type) != 0
                    and allele_freq >= self.params.basic_filter.min_mut_freq):
                    
                    # define allele
                    Abase, Tbase, Gbase, Cbase = self.retrieve_base_string_each_base_type(A,T,G,C)
                    _all_base = Abase + Gbase + Cbase + Tbase
                    alt = AlignmentReadsStats.define_allele(_all_base, ref=ref_base)
                
                    # Array in seq with read strand information
                    #reads_strand_information(G, 'rev')
                    G_base_r = [_.alignment.seq[_.qpos] for _ in G
                                if _.alignment.is_reverse]
                    G_base_f = [_.alignment.seq[_.qpos] for _ in G
                                if not _.alignment.is_reverse]
                    
                    A_base_r = [_.alignment.seq[_.qpos] for _ in A
                                if _.alignment.is_reverse]
                    A_base_f = [_.alignment.seq[_.qpos] for _ in A
                                if not _.alignment.is_reverse]
                
                    T_base_r = [_.alignment.seq[_.qpos] for _ in T
                                if _.alignment.is_reverse]
                    T_base_f = [_.alignment.seq[_.qpos] for _ in T
                                if not _.alignment.is_reverse]
                    
                    C_base_r = [_.alignment.seq[_.qpos] for _ in C
                                if _.alignment.is_reverse]
                    C_base_f = [_.alignment.seq[_.qpos] for _ in C
                                if not _.alignment.is_reverse]

                    dp4 = (AlignmentReadsStats.compute_dp4(ref_base,
                                            len(A_base_r), len(A_base_f),
                                            len(T_base_r), len(T_base_f),
                                            len(G_base_r), len(G_base_f),
                                            len(C_base_r), len(C_base_f)))
                    ## TODO:
                    ## comapre speed by __len__() and count()
                    #Ac = [_.alignment.seq[_.qpos] for _ in A].count('A')
                    #Tc = [_.alignment.seq[_.qpos] for _ in T].count('T')
                    #Gc = [_.alignment.seq[_.qpos] for _ in G].count('G')
                    # G_base_count = G_base_r + G_base_f
                    #Cc = [_.alignment.seq[_.qpos] for _ in C].count('C')
                    
                    ###########################
                    ### Statistical filsher ###
                    ###########################
                    strand_bias_p = float()
                    positional_bias_p = float()
                    base_call_bias_p = float()
                    if (self.params.stat_filter.strand_bias):
                        strand_bias_p = strand_bias_filter(passed_matches, passed_mismatches)
                        if strand_bias_p > self.params.stat_filter.sig_level:
                            pass
                            
                    if (self.params.stat_filter.pos_bias):
                        positional_bias_p = 0
                    if (self.params.stat_filter.baq_bias):
                        base_call_bias_p = 0
                        
                    #if pos == 47721228:
                    #    print ref_base, pos
                    #    print  coverage
                    #    #qpos = [_.alignment.qual[_.qpos] for _ in passed_reads]
                    # 
                    #    mis_r = [_.alignment.seq[_.qpos] for _ in passed_mismatches if _.alignment.is_reverse]
                    #    mis_f = [_.alignment.seq[_.qpos] for _ in passed_mismatches if not _.alignment.is_reverse]
                    #    ma_r =  [_.alignment.seq[_.qpos] for _ in passed_matches if _.alignment.is_reverse]
                    #    ma_f =  [_.alignment.seq[_.qpos] for _ in passed_matches if not _.alignment.is_reverse]
                    #    
                    #    print [_.qpos for _ in passed_reads]
                    #    print [_.alignment.qend for _ in passed_mismatches]
                    #    print [_.alignment.qstart for _ in passed_reads]
                    #    #raise SystemExit
                    #    
                    #if (positional_bias_p > self.params.stat_filter.sig_level
                    #    or base_call_bias_p > self.params.stat_filter.sig_level
                    #    or strand_bias_p > self.params.stat_filter.sig_level):
                    if 1: 
                        yield {
                            'chrom': bam_chrom,
                            'pos': pos,
                            'ref': ref_base,
                            'alt': alt,
                            'coverage': len(passed_reads),
                            'mismatches': len(passed_mismatches),
                            'matches': len(passed_matches),
                            'cov': coverage,
                            'mismatch_ratio': allele_freq,
                            'types': mutation_type,
                            'dp4': dp4
                        }
                        
    def retrieve_reads_with_specify_base(self, reads, base):
        '''
        Args:
         pysam.csamtools.PileupRead(object)
         Specify nucleotide base type(str)
        Returns:
         pysam.csamtools.PileupRead object
        '''
        
        if base not in  ['A', 'T', 'G', 'C']:
            raise ValueError()
        return [_ for _ in reads if _.alignment.seq[_.qpos] == base]

    def retrieve_reads_each_base_type(self, reads):
        '''
        Args:
         pysam.csamtools.PileupRead objects(list)
        Returns:
         pysam.csamtools.PileupRead object
        '''
        
        A = [_ for _ in reads if _.alignment.seq[_.qpos] == 'A']
        T = [_ for _ in reads if _.alignment.seq[_.qpos] == 'T']
        C = [_ for _ in reads if _.alignment.seq[_.qpos] == 'C']
        G = [_ for _ in reads if _.alignment.seq[_.qpos] == 'G']
        return A, T, G, C

    def retrieve_base_string_each_base_type(self, A, T, G, C):
        '''
        Array in four type of sequence
        
        Args:
         pysam.csamtools.PileupRead object(list)
        Returns:
         4type of base strings(list)
        '''
        Gb =  [_.alignment.seq[_.qpos] for _ in G]
        Ab =  [_.alignment.seq[_.qpos] for _ in A]
        Tb =  [_.alignment.seq[_.qpos] for _ in T]
        Cb =  [_.alignment.seq[_.qpos] for _ in C]
        return Ab, Gb, Cb, Tb

    def retrieve_base_string_with_strand(self, base, strand=None):
        '''
        Args:
         base(list), e.g. data[pysam.csamtools.PileupRead object]
         strand=1 or 0. 1 is forward strand, 0 is reverse strand.
        
        Returns:
         bases with specify strand
        '''
        
        forward = [_.alignment.seq[_.qpos] for _ in G
                   if _.alignment.is_reverse]
        reverse = [_.alignment.seq[_.qpos] for _ in G
                   if not _.alignment.is_reverse]
        if strnad == 1:
            return forward
        elif strand == 0:
            return reverse
        else:
            return None
            
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

    
    def fasta_info(self):
        # info. for fasta
        print "### info. for fasfile object ###"
        print "filename: %s" % self.fafile.filename

    def sam_info(self):
        # info. for loaded samfile
        print "### info. for samfile object from given Bam header @SQ ###"
        print "Sam file: %s" % self.samfile.filename
        print "lengths: %s" % [_ for _ in self.samfile.lengths]
        print "mapped %d: " % self.samfile.mapped
        print "N_references: %s" % self.samfile.nreferences
        print "references: %s" % [_ for _ in self.samfile.references]
        print "unmapped: %s" % self.samfile.unmapped


def _parse_faidx(filename):
    fasta_chrom_name = []
    with open(filename, 'r') as fh:
        for row in fh:
            data = row.split('\t')
            fasta_chrom_name.append(data[0])
    return fasta_chrom_name
            
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
            raise RuntimeError('{filename:s} of faidx file is not found'.
                               format(filename=fa_dx_filename))
            
def _resolve_chrom_name(bam_chr=None, fa_chr=None):
    raise NotImplementedError()
    if not fa_chr.startswith('chr'):
        return 'chr' + fa_chr
    else:
        return fa_chr
            
if __name__ == '__main__':
    align = AlignmentStream("params")
    
    
    
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
    
class AlignmentStreamMerger(object):
    def __init__(self, rna, dna):
        raise NotImplementedError()
        self.rna = rna
        self.dna = dna

    def merge_streaming(self):
        dna_stream = AlignmentStream(conf)
        rna_stream = AlignmentStream(conf)

