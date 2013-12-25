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
from Ivy.alignment.filters import strand_bias_filter, positional_bias_filter
from Ivy.alignment.stats import AlignmentReadsStats

__program__ = 'stream'
__author__ = 'Soh Ishiguro <yukke@g-language.org>'
__license__ = ''
__status__ = 'development'

class BaseStringGenerator(object):
    @classmethod
    def retrieve_base_string_each_base_type(self, a=[], t=[], g=[], c=[]):
        '''
        Array in four type of sequence
        
        Args:
         pysam.csamtools.PileupRead object(list)
        Returns:
         4type of base strings(dict)
        '''
        
        Gb =  [_.alignment.seq[_.qpos] for _ in g]
        Ab =  [_.alignment.seq[_.qpos] for _ in a]
        Tb =  [_.alignment.seq[_.qpos] for _ in t]
        Cb =  [_.alignment.seq[_.qpos] for _ in c]
        return {'A': Ab, 'T': Tb, 'G': Gb, 'C': Cb}
            
    @classmethod
    def retrieve_base_string_with_strand(self, base, strand=None):
        '''
        Args:
         base(list), e.g. data[pysam.csamtools.PileupRead object]
         strand=1 or 0. 1 is forward strand, 0 is reverse strand.
        
        Returns:
         bases with specify strand
        '''
        
        if strand == 1:
            return [_.alignment.seq[_.qpos] for _ in base
                    if _.alignment.is_reverse]
        elif strand == 0:
            return [_.alignment.seq[_.qpos] for _ in base
                    if not _.alignment.is_reverse]
        else:
            return None

            
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
                      and not _.alignment.is_secondary
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

DEBUG = False
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

        ### Resolve to explore specific region or not
        # Flag == 1: explore all region
        if self.params.region.all_flag == 1:
            self.params.region.start = None
            self.params.region.end = None
            self.params.region.chrom = None

        # Flag == 0: explore specific region
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
            ### Loading alingnment with params ###
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
                    raise SystemExit("Called methods: '{0:s}'".format(self.reads_allow_duplication.__name__))
                passed_reads, passed_matches, passed_mismatches = self.reads_allow_duplication(col.pileup, ref_base)
                
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
                    raise SystemExit("Called methods: '{0:s}'".format(self.reads_without_filter.__name__))
                passed_reads, passed_matches, passed_mismatches = self.reads_without_filter(col.pileups, ref_base)
            
            ##############################
            ### Basic filters in reads ###
            ##############################
            # --min-rna-baq
            alignstat = AlignmentReadsStats()
            quals_in_pos = alignstat.quals_in_pos(passed_reads)
            average_baq = alignstat.average_base_quality(passed_reads)
            
            # --min-rna-cov
            coverage = alignstat.reads_coverage(passed_reads)
            
            # --min-rna-mapq
            average_mapq = alignstat.average_mapq(passed_reads)
            
            specific_reads = self.retrieve_reads_each_base_type(passed_reads)
            A_reads = specific_reads.get('A')
            T_reads = specific_reads.get('T')
            G_reads = specific_reads.get('G')
            C_reads = specific_reads.get('C')
            
            if (self.params.basic_filter.min_rna_cov <= coverage
                and self.params.basic_filter.min_rna_mapq <= average_mapq
                and self.params.basic_filter.min_baq_rna <= average_baq):
                
                # --min-mis-frequency
                allele_freq = alignstat.mismatch_frequency(m=passed_matches, mis=passed_mismatches)
                ag_freq = alignstat.a_to_g_frequency(a=A_reads, g=G_reads)
                
                # --num-allow-type
                mutation_type = self.mutation_types(A_reads, T_reads, G_reads, C_reads, ref=ref_base)
                
                if (len(mutation_type) <= self.params.basic_filter.num_type
                    and len(mutation_type) != 0
                    and allele_freq >= self.params.basic_filter.min_mut_freq):
                    
                    basegen = BaseStringGenerator()
                    base = basegen.retrieve_base_string_each_base_type(a=A_reads, t=T_reads, g=G_reads, c=C_reads)
                    
                    Abase = base.get('A')
                    Tbase = base.get('T')
                    Gbase = base.get('G')
                    Cbase = base.get('C')
                    
                    _all_base = Abase + Gbase + Cbase + Tbase
                    alt = alignstat.define_allele(_all_base, ref=ref_base)
                    
                    # Specific base string by read strand(forward/reverse)
                    G_base_r = basegen.retrieve_base_string_with_strand(G_reads, strand=0)
                    G_base_f = basegen.retrieve_base_string_with_strand(G_reads, strand=1)
                    
                    A_base_r = basegen.retrieve_base_string_with_strand(A_reads, strand=0)
                    A_base_f = basegen.retrieve_base_string_with_strand(A_reads, strand=1)
                    
                    T_base_r = basegen.retrieve_base_string_with_strand(T_reads, strand=0)
                    T_base_f = basegen.retrieve_base_string_with_strand(T_reads, strand=1)
                    
                    C_base_r = basegen.retrieve_base_string_with_strand(C_reads, strand=0)
                    C_base_f = basegen.retrieve_base_string_with_strand(C_reads, strand=1)
                    
                    dp4 = (alignstat.compute_dp4(ref=ref_base,
                                                 ar=len(A_base_r), af=len(A_base_f),
                                                 tr=len(T_base_r), tf=len(T_base_f),
                                                 gr=len(G_base_r), gf=len(G_base_f),
                                                 cr=len(C_base_r), cf=len(C_base_f)))
                    
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

                    positional_bias_p = positional_bias_filter(m=passed_matches, mis=passed_mismatches)
                    if (self.params.stat_filter.pos_bias):
                        pass
                        
                    if (self.params.stat_filter.baq_bias):
                        base_call_bias_p = 0
                        
                    # Faital error if diff. in len(N) != (len(Nr)+len(Nf))
                    # TODO: Wrapp *Error class in error.py
                    if len(Abase) != len(A_base_r + A_base_f):
                        raise ValueError, ("All{all:0}, for{f:1}, rev{r:1}".format(
                            all=len(Abase), f=len(A_base_f), r=len(A_base_r)))
                        
                    if len(Tbase) != len(T_base_r + T_base_f):
                        raise ValueError, ("All{all:0}, for{f:1}, rev{r:1}".format(
                            all=len(Tbase), f=len(T_base_f), r=len(T_base_r)))
                        
                    if len(Gbase) != len(G_base_r + G_base_f):
                        raise ValueError, ("All{all:0}, for{f:1}, rev{r:1}".format(
                            all=len(Gbase), f=len(G_base_f), r=len(G_base_r)))

                    if len(Cbase) != len(C_base_r + C_base_f):
                        raise ValueError, ("All: {all:0}, for: {f:1}, rev: {r:1} at {pos:2}, {ref:3}, {m:4}".format(
                            all=len(Cbase), f=len(C_base_f), r=len(C_base_r), pos=pos, ref=ref_base, m=mutation_type))

                    #assert False, ([_.alignment.tlen for _ in passed_mismatches])
                    
                    d = {
                        'chrom': bam_chrom,
                        'pos': pos,
                        'ref': ref_base,
                        'alt': alt[0],
                        'coverage': len(passed_reads),
                        'mismatches': len(passed_mismatches),
                        'matches': len(passed_matches),
                        'allele_freq': allele_freq,
                        'positiona_bias': positional_bias_p,
                        #'ag_freq': ag_freq,
                        #'types': mutation_type,
                        #'dp4': dp4,
                        #'average_baq': average_baq,
                        #'average_mapq': average_mapq,
                        #qual_in_pos': quals_in_pos,
                        #'raw_quals': [_.alignment.qual[_.qpos] for _ in passed_reads],
                        #'mutation_type': mutation_type,
                        #'A': Abase,
                        #'G': Gbase,
                        #'T': Tbase,
                        #'C': Cbase,
                        #'A_f': A_base_f,
                        #'A_r': A_base_r,
                        #'G_f': G_base_f,
                        #'G_r': G_base_r,
                        #'T_f': T_base_f,
                        #'T_r': T_base_r,
                        #'C_f': C_base_f,
                        #'C_r': C_base_r,
                        }
                    yield d

                        
    def mutation_types(self, A, T, G, C, ref=None):
        mutation_types = {}
        if len(A) > 0 and ref != 'A':
            mutation_types.update({'A': len(A)})
        elif len(T) > 0 and ref != 'T':
            mutation_types.update({'T': len(T)})
        elif len(G) > 0 and ref != 'G':
            mutation_types.update({'G': len(G)})
        elif len(C) > 0 and ref != 'C':
            mutation_types.update({'C': len(C)})
        return mutation_types
        
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
        
        a = [_ for _ in reads if _.alignment.seq[_.qpos] == 'A']
        t = [_ for _ in reads if _.alignment.seq[_.qpos] == 'T']
        c = [_ for _ in reads if _.alignment.seq[_.qpos] == 'C']
        g = [_ for _ in reads if _.alignment.seq[_.qpos] == 'G']
        return {'A': a, 'T': t, 'G': g, 'C': c}
            
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
            
    def _resolve_chrom_name(self, bam_chr=None, fa_chr=None):
        raise NotImplementedError()
        if not fa_chr.startswith('chr'):
            return 'chr' + fa_chr
        else:
            return fa_chr

        
if __name__ == '__main__':
    align = AlignmentStream("params")


