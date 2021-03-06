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
from Ivy.utils import die, AttrDict, IvyLogger, convert_base
from Ivy.alignment.filters import strand_bias_filter, positional_bias_filter
from Ivy.alignment.stats import AlignmentReadsStats
from Ivy.annotation.annotation import GTF

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
            # reverse
            return [_.alignment.seq[_.qpos] for _ in base
                    if _.alignment.is_reverse]
        elif strand == 0:
            # forward
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

    def reads_filter_by_all_params(self, pileup_obj, ref_base, is_single=False):
        '''
        List of read filters:
         Proper pair reads
         NOT duplicated reads
         NOT unmapped reads
         NOT deletion reads, note that is_del is quite effective for reducing reads
         NOT fail to quality check reads
         NOT paased spliced reads (RNA-seq data)
         NOT passed primary alignmnet reads
        '''
        
        _reads = []
        # paired-end read (default)
        if not is_single:
            _reads = [_ for _ in pileup_obj
                      if (_.alignment.is_proper_pair
                          and not _.alignment.is_qcfail
                          and not _.alignment.is_duplicate
                          and not _.alignment.is_unmapped
                          and not _.alignment.is_secondary
                          and _.indel == 0
                          and _.is_del == 0)]
        # single-end
        else:
            _reads = [_ for _ in pileup_obj
                      if (not _.alignment.is_qcfail
                          and not _.alignment.is_duplicate
                          and not _.alignment.is_unmapped
                          and not _.alignment.is_secondary
                          and _.indel == 0
                          and _.is_del == 0)]
            
        _matches = [_ for _ in _reads if _.alignment.seq[_.qpos] == ref_base]
        _mismatches = [_ for _ in _reads if _.alignment.seq[_.qpos] != ref_base]
        return _reads, _matches, _mismatches
            
    def reads_filter_without_pp(self, pileup_obj, ref_base):
        '''
        List of read filters:
         NOT fail to quality check reads
         NOT duplicated reads
         NOT unmapped reads
         NOT deletion reads
        '''
        
        _reads = [_ for _ in pileup_obj
                  if (_.alignment.is_qcfail
                      and not _.alignment.is_duplicate
                      and not _.alignment.is_unmapped
                      and _.is_del == 0
                      and _.indel == 0)]
        
        _matches = [_ for _ in _reads if _.alignment.seq[_.qpos] == ref_base]
        _mismatches = [_ for _ in _reads if _.alignment.seq[_.qpos] != ref_base]
        return _reads, _matches, _mismatches
        
    def reads_allow_duplication(self, pileup_obj, ref_base):
        '''
        List of read filters:
         proper pair reads
         NOT failt to quality check reads
         NOT unmapped reads
         NOT deletion reads
        '''
        
        _reads = [_ for _ in col.pileup_obj
                  if (_.alignment.is_proper_pair
                      and not _.alignment.is_qcfail
                      and not _.alignment.is_unmapped
                      and not _.is_del
                      and not _.indel == 0)]
        
        _matches = [_ for _ in _reads if _.alignment.seq[_.qpos] == ref_base]
        _mismatches = [_ for _ in _reads if _.alignment.seq[_.qpos] != ref_base]
        return _reads, _matches, _mismatches
            
    def reads_without_filter(self, pileup_obj, ref_base):
        '''All reads are passed through under this filter'''
        
        _reads = [_ for _ in pileup_obj
                  if (_.indel == 0
                      and _.is_del == 0)]
        
        _matches = [_ for _ in _reads if _.alignment.seq[_.qpos] == ref_base]
        _mismatches = [_ for _ in _reads if _.alignment.seq[_.qpos] != ref_base]
        return _reads, _matches, _mismatches
        
    def reads_allow_insertion(self):
        print "Use --rm-insertion-reads is recommended"
        return None
        
    def reads_allow_deletion(self, pileup_obj, ref_base):
        print "Use --rm-deletion-reads is recommended"
        return None

DEBUG = False
reads_filter_params_debug = False
class AlignmentStream(FilteredAlignmentReadsGenerator):
    def __init__(self, __params, mol_type=None):
        '''
        Initialize for pileup bam files to explore RDD sites
        
        Args:
         Alignment parameters wrapped by AttrDic object
         mol_type(str)='RNA' or 'DNA'
        
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

            
        if mol_type == 'RNA':
            self.samfile = self.__load_bam(rna=True)
        elif mol_type == 'DNA':
            self.samfile = self.__load_bam(rna=False)
        else:
            raise ValueError("Molecular type(DNA/RNA) is None! then you specify mol_type=''")
        
        self.fafile = pysam.Fastafile(self.params.fasta)
        self.one_based = self.params.one_based
        
        ### Resolve to explore specific region or not
        # explore specific chrom each thred if multi threading
        if self.params.n_threads:
            self.params.region.all_flag = 0
            self.params.region.start = None
            self.params.region.end = None
            
        # Flag == 1: explore all region
        if self.params.region.all_flag == 1:
            self.params.region.start = None
            self.params.region.end = None
            self.params.region.chrom = None
            
        # Flag == 0: explore specific region
        elif self.params.region.all_flag == 0:
            if (self.params.region.chrom
                and self.params.region.start and self.params.region.end):
                self.params.region.start, self.params.region.end = \
                                                                   self.__resolve_coords(
                                                                       self.params.region.start,
                                                                       self.params.region.end,
                                                                       self.params.one_based)
        else:
            # explore all region if self.params.region.* is None
            pass
                
        #if self.params.region.chrom.startswith('chr'):
        #    self.params.region.chrom = self.params.region.chrom.replace('chr', '', 1)
        #else:
        #    self.params.region.chrom = self.params.region.chrom

        if self._is_same_chromosome_name(bam=self.params.r_bams, fa=self.params.fasta):
            pass
        else:
            #raise RuntimeError("Invalid chromosome name in {fa}, {reg}".format(fa=self.params.fasta, reg=self.params.region))
            return None
            
        if self.params.verbose:
            self.logger.debug(AttrDict.show(self.params))
            
        if DEBUG:
            _fasta_info()
            _sam_info()

        self.gtf_cls = GTF(self.params.gtf)
        

    def __load_bam(self, rna=True):
        if rna:
            return pysam.Samfile(self.params.r_bams, 'rb', check_header=True, check_sq=True)
        elif not rna:
            return pysam.Samfile(self.params.d_bams, 'rb', check_header=True, check_sq=True)
            
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
        '''
        Returns:
         Pileuped alignment data(tuple):
          passsed_reads[0],
          passed_matches[1],
          passed_mismatches[2]

        '''
        
        if self.params.verbose:
            self.logger.debug("Start pileup bam '{0}' file...".format(self.__class__.__name__))

        try:
            for col in self.samfile.pileup(reference=self.params.region.chrom,
                                           start=self.params.region.start,
                                           end=self.params.region.end):
                self.bam_chrom = self.samfile.getrname(col.tid)
                
                if self.params.one_based:
                    self.pos = col.pos + 1
                else:
                    self.pos = col.pos
                self.ref_base= self.fafile.fetch(reference=self.bam_chrom,
                                            start=col.pos,
                                            end=col.pos+1).upper()
                if not self.ref_base:
                    # Skip if chrom in bam is not in reference genome
                    continue
                    #raise ValueError(
                    #   'No sequence content within chroms: {chrom:s}, start: {start:s}, end: {end:s}'.format(
                    #       chrom=self.params.region.chrom, start=self.params.region.start, end=self.params.region.end))
                elif self.ref_base == 'N' or self.ref_base == 'n':
                    continue
                 
                # filter reads with all params
                
                if (self.params.basic_filter.rm_duplicated
                    and self.params.basic_filter.rm_deletion
                    and self.params.basic_filter.rm_insertion):
                    if reads_filter_params_debug:
                        self.params.show(self.params.basic_filter)
                        raise SystemExit("Method: {0:s}".format(
                            self.reads_filter_by_all_params.__name__))
                        
                    passed_reads, passed_ma, passed_mis = (
                        self.reads_filter_by_all_params(col.pileups, self.ref_base, is_single=self.params.is_single))
                    yield {'reads': passed_reads,
                           'ma': passed_ma,
                           'mis': passed_mis}
                    
                # allow duplicated containing reads
                elif (not self.params.basic_filter.rm_duplicated
                      and self.params.basic_filter.rm_deletion
                      and self.params.basic_filter.rm_insertion):
                    if reads_filter_params_debug:
                        self.params.show(self.params.basic_filter)
                        raise SystemExit("Method: '{0:s}'".format(
                            self.reads_allow_duplication.__name__))
                        
                    passed_reads, passed_ma, passed_mis = (
                        self.reads_allow_duplication(col.pileups, self.ref_base))
                    yield {'reads': passed_reads,
                           'ma': passed_ma,
                           'mis': passed_mis}
                    
                # allow deletions containing reads
                elif (not self.params.basic_filter.rm_deletion
                      and self.params.basic_filter.rm_insertion
                      and self.params.basic_filter.rm_duplicated):
                    if reads_filter_params_debug:
                        self.params.show(self.params.basic_filter)
                        raise SystemExit("Method: {0:s}".format(self.reads_allow_deletion.__name__))
                        
                    passed_reads, passed_ma, passed_mis = (
                        self.reads_allow_deletion(col.pileups, self.ref_base))
                    yield {'reads': passed_reads,
                           'ma': passed_ma,
                           'mis': passed_mis}
                    
                # allow insertion containing reads
                elif (not self.params.basic_filter.rm_insertion
                      and self.params.basic_filter.rm_deletion
                      and self.params.basic_filter.rm_duplicated):
                    if reads_filter_params_debug:
                        self.params.show(self.params.basic_filter)
                        raise SystemExit("Method: {0:s}".format(
                            self.params.basic_filter, self.reads_allow_insertion.__name__))
                        
                    passed_reads, passed_ma, passed_mis = (
                        self.reads_allow_insertion(col.pileups, self.ref_base))
                    yield {'reads': passed_reads,
                           'ma': passed_ma,
                           'mis': passed_mis}
     
                # no filter
                else:
                    if reads_filter_params_debug:
                        self.params.show(self.params.basic_filter)
                        raise SystemExit("Method: '{0:s}'".format(self.reads_without_filter.__name__))
                        
                    passed_reads, passed_ma, passed_mis = (
                        self.reads_without_filter(col.pileups, self.ref_base))
                    yield {'reads': passed_reads,
                           'ma': passed_ma,
                           'mis': passed_mis}
                
        except ValueError:
            # e.g. ValueError: invalid reference `chr18` in pysam error
            # return tuple in None
            #yield tuple([None, None, None])
            yield {'reads': None, 'ma': None, 'mis': None}
            
    def mutation_types(self, A=None, T=None, G=None, C=None, ref=None):
        '''
        Define possible mutation types in given list in each base type
        '''
        mutation_types = {}
        if len(A) > 0 and ref != 'A':
            mutation_types.update({'A': len(A)})
        if len(T) > 0 and ref != 'T':
            mutation_types.update({'T': len(T)})
        if len(G) > 0 and ref != 'G':
            mutation_types.update({'G': len(G)})
        if len(C) > 0 and ref != 'C':
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
        '''info. for fasta'''
        
        print "### info. for fasfile object ###"
        print "filename: %s" % self.fafile.filename

    def sam_info(self):
        '''info. for loaded samfile'''
        
        print "### info. for samfile object from given Bam header @SQ ###"
        print "Sam file: %s" % self.samfile.filename
        print "lengths: %s" % [_ for _ in self.samfile.lengths]
        print "mapped %d: " % self.samfile.mapped
        print "N_references: %s" % self.samfile.nreferences
        print "references: %s" % [_ for _ in self.samfile.references]
        print "unmapped: %s" % self.samfile.unmapped


    def _parse_faidx(self, filename):
        '''
        Return:
         chromosomes(list) in given fai file
        
        Example:
         $ cat genome.fa.fai
         chr18 78077248 7 50 51
         chr19 59128983 79638807 50 51
         chr20 63025520 139950377 50 51
        '''
        
        fasta_chrom_name = []
        with open(filename, 'r') as fh:
            for row in fh:
                data = row.split('\t')
                fasta_chrom_name.append(data[0])
        return fasta_chrom_name
                
    def _is_same_chromosome_name(self, bam=None, fa=None):
        __bam = pysam.Samfile(os.path.abspath(bam), 'rb')
        __fa = pysam.Fastafile(os.path.abspath(fa))
        bam_references = __bam.references
        fa_filename = __fa.filename
        fa_dx_filename = fa_filename + '.fai'
        
        if os.path.isfile(fa_dx_filename):
            for bam_chr in bam_references:
                if any([fai_chr == bam_chr for fai_chr in self._parse_faidx(fa_dx_filename)]):
                    return True
                else:
                    #print tuple([self._parse_faidx(fa_dx_filename), bam_references])
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
            
class RNASeqAlignmentStream(AlignmentStream):
    def __init__(self, rna_params):
        AlignmentStream.__init__(self, rna_params, mol_type='RNA')
        self.params = rna_params
    
    def filter_stream(self):
        for data in self.pileup_stream():
            passed_mismatches = data['mis']
            if passed_mismatches < 2:
                continue
            passed_reads = data['reads']
            passed_matches = data['ma']

            debug_pos = False
            
            if self.pos == debug_pos:
                specific_reads = self.retrieve_reads_each_base_type(passed_reads)
                A_reads = specific_reads.get('A')
                T_reads = specific_reads.get('T')
                G_reads = specific_reads.get('G')
                C_reads = specific_reads.get('C')
                base = basegen.retrieve_base_string_each_base_type(a=A_reads,
                                                                   t=T_reads,
                                                                   g=G_reads,
                                                                   c=C_reads)
                Abase = base.get('A')
                Tbase = base.get('T')
                Gbase = base.get('G')
                Cbase = base.get('C')
                print "A: ", Abase
                print "T: ", Tbase
                print "G: ", Gbase
                print "C: ", Cbase
                
            ##############################
            ### Basic filters in reads ###
            ##############################
            alignstat = AlignmentReadsStats()
            
            # --min-rna-cov
            coverage = alignstat.reads_coverage(passed_reads)
            if coverage <= self.params.basic_filter.min_rna_cov:
                continue
                
            # --min-rna-baq                
            average_baq = alignstat.average_base_quality(passed_reads)
            baq_in_mismatch = alignstat.average_base_quality(passed_mismatches)
            if baq_in_mismatch <= self.params.basic_filter.min_rna_baq:
                continue
                
            if self.pos == debug_pos:
                print "base qual: {0}, baq in mis: {1}".format(average_baq, baq_in_mismatch)
            
            # --min-rna-mapq
            average_mapq = alignstat.average_mapq(passed_reads)
            if self.pos == debug_pos:
                print "mapq: {0}".format(average_mapq)
            if average_mapq <= self.params.basic_filter.min_rna_mapq:
                continue

            # --min-mis-frequency
            allele_freq = alignstat.mismatch_frequency(m=passed_matches, mis=passed_mismatches)
            if self.pos == debug_pos:
                print "Allele freq. {0}".format(allele_freq)
            if allele_freq <= self.params.basic_filter.min_mut_freq:
                continue
                
            # --num-allow-type
            specific_reads = self.retrieve_reads_each_base_type(passed_reads)
            A_reads = specific_reads.get('A')
            T_reads = specific_reads.get('T')
            G_reads = specific_reads.get('G')
            C_reads = specific_reads.get('C')

            mutation_type = self.mutation_types(A=A_reads, T=T_reads, G=G_reads, C=C_reads, ref=self.ref_base)
            if self.pos == debug_pos:
                print "Mut. type: {0}".format(mutation_type)
            
            if len(mutation_type) == 0:
                continue
                
            if len(mutation_type) > self.params.basic_filter.num_type:
                continue

            basegen = BaseStringGenerator()
            base = basegen.retrieve_base_string_each_base_type(a=A_reads,
                                                               t=T_reads,
                                                               g=G_reads,
                                                               c=C_reads)
            Abase = base.get('A')
            Tbase = base.get('T')
            Gbase = base.get('G')
            Cbase = base.get('C')
                    
            _all_base = Abase + Gbase + Cbase + Tbase
            _alt = alignstat.define_allele(_all_base, ref=self.ref_base)
            alt_base = ",".join([",".join(_[0]) for _ in _alt])
            
            # Specific base string by read strand(forward/reverse)
            G_base_r = basegen.retrieve_base_string_with_strand(G_reads, strand=0)
            G_base_f = basegen.retrieve_base_string_with_strand(G_reads, strand=1)
                    
            A_base_r = basegen.retrieve_base_string_with_strand(A_reads, strand=0)
            A_base_f = basegen.retrieve_base_string_with_strand(A_reads, strand=1)
                    
            T_base_r = basegen.retrieve_base_string_with_strand(T_reads, strand=0)
            T_base_f = basegen.retrieve_base_string_with_strand(T_reads, strand=1)
                    
            C_base_r = basegen.retrieve_base_string_with_strand(C_reads, strand=0)
            C_base_f = basegen.retrieve_base_string_with_strand(C_reads, strand=1)

            ###########################
            ### Statistical filsher ###
            ###########################
            # and False means p<sig_level location
            strand_bias_p = .0
            if self.params.stat_filter.strand_bias:
                strand_bias_p = strand_bias_filter(m=passed_matches,
                                                   mis=passed_mismatches)
                if strand_bias_p < self.params.stat_filter.sig_level:
                    continue
            positional_bias_p = .0
            if (self.params.stat_filter.pos_bias):
                positional_bias_p = positional_bias_filter(m=passed_matches,
                                                                   mis=passed_mismatches)
                if positional_bias_p < self.params.stat_filter.sig_level:
                    continue
            base_call_bias_p = .0
            if (self.params.stat_filter.baq_bias):
                base_call_bias_p = .0
                if base_call_bias_p < self.params.stat_filter.sig_level:
                    continue
                    
            dp4 = (alignstat.compute_dp4(ref=self.ref_base,
                                         ar=len(A_base_r), af=len(A_base_f),
                                         tr=len(T_base_r), tf=len(T_base_f),
                                         gr=len(G_base_r), gf=len(G_base_f),
                                         cr=len(C_base_r), cf=len(C_base_f)))
            ag_freq = alignstat.a_to_g_frequency(a=Abase, g=Gbase)

            # Considering expressed strand in transcript
            _strand = self.gtf_cls.strand_info(self.bam_chrom, start=self.pos)
            if _strand is not None or _strand == "Intergenic":
                #print "[Original] %s-%s, ref: %s, alt: %s, strand: %s" % (self.bam_chrom, self.pos, self.ref_base, alt_base, _strand),
                _type = convert_base(ref=self.ref_base, alt=alt_base, strand=_strand)
                _ref = _type[0]
                _alt = _type[1]
                alt_base = _alt
                #print "[Updated] ref: %s, alt: %s, strand: %s" % (_ref, _alt, _strand)
            else:
                _type = self.ref_base + "-to-" + alt_base
            
            d = {
                'chrom': self.bam_chrom,
                'pos': self.pos,
                'id': '.',
                'ref': self.ref_base,
                'alt': alt_base,
                'qual': baq_in_mismatch,
                'filt': '.',
                'coverage': len(passed_reads),
                'mismatches': len(passed_mismatches),
                'matches': len(passed_matches),
                'allele_freq': allele_freq,
                'positional_bias': positional_bias_p,
                'strand_bias': strand_bias_p,
                'type': _type,
                'base_call_bias': base_call_bias_p,
                'ag_freq': ag_freq,
                'types': mutation_type,
                'dp4': dp4,
                'average_baq': average_baq,
                'average_mapq': average_mapq,
                'mutation_type': mutation_type,
                'A': Abase,
                'G': Gbase,
                'T': Tbase,
                'C': Cbase,
                'A_f': A_base_f,
                'A_r': A_base_r,
                'G_f': G_base_f,
                'G_r': G_base_r,
                'T_f': T_base_f,
                'T_r': T_base_r,
                'C_f': C_base_f,
                'C_r': C_base_r,

            }
            yield d
            
            
                        
class DNASeqAlignmentStream(AlignmentStream):
    def __init__(self, dna_params):
        AlignmentStream.__init__(self, dna_params, mol_type='DNA')
        self.params = dna_params
        
    def filter_stream(self):
        for data in self.pileup_stream():
            passed_mismatches = data[2]
            if len(passed_mismatches) < 1:
                continue
            passed_reads = data[0]
            passed_matches = data[1]
            
            alignstat = AlignmentReadsStats()
            
            # --min-dna-cov
            coverage = alignstat.reads_coverage(passed_reads)
            if coverage <= self.params.basic_filter.min_dna_cov:
                continue
                
            # --min-dna-baq
            average_baq = alignstat.average_base_quality(passed_reads)
            if average_baq <= self.params.basic_filter.min_dna_baq:
                continue
                
            # --min-dna-mapq
            average_mapq = alignstat.average_mapq(passed_reads)
            if average_mapq <= self.params.basic_filter.min_dna_mapq:
                continue
                
            # --min-mis-frequency (it just mean as allele frequency in SNPs)
            allele_freq = alignstat.mismatch_frequency(m=passed_matches, mis=passed_mismatches)
            if allele_freq <= self.params.basic_filter.min_mut_freq:
                continue

            # --num-allow-type
            specific_reads = self.retrieve_reads_each_base_type(passed_reads)
            A_reads = specific_reads.get('A')
            T_reads = specific_reads.get('T')
            G_reads = specific_reads.get('G')
            C_reads = specific_reads.get('C')
            mutation_type = self.mutation_types(A=A_reads, T=T_reads, G=G_reads, C=C_reads, ref=self.ref_base)
            if len(mutation_type) == 0:
                continue
            
            basegen = BaseStringGenerator()
            base = basegen.retrieve_base_string_each_base_type(a=A_reads,
                                                               t=T_reads,
                                                               g=G_reads,
                                                               c=C_reads)
            Abase = base.get('A')
            Tbase = base.get('T')
            Gbase = base.get('G')
            Cbase = base.get('C')
            _all_base = Abase + Gbase + Cbase + Tbase
            _alt = alignstat.define_allele(_all_base, ref=self.ref_base)
            alt_base = ",".join([",".join(_[0]) for _ in _alt])
                
            d = {
                'chrom': self.bam_chrom,
                'pos': self.pos,
                'ref': self.ref_base,
                'alt': alt_base,
                #'coverage': len(passed_reads),
                #'mismatches': len(passed_mismatches),
                #'matches': len(passed_matches),
                #'allele_freq': allele_freq,
                #'positional_bias': positional_bias_p,
                #'strand_bias': strand_bias_p,
                #'base_call_bias': base_call_bias_p,
                #'ag_freq': ag_freq,
                #'types': mutation_type,
                #'dp4': dp4,
                #'average_baq': average_baq,
                #'average_mapq': average_mapq,
                #'qual_in_pos': quals_in_pos,
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
