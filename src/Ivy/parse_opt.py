from optparse import OptionParser, OptionGroup, HelpFormatter, IndentedHelpFormatter
import os.path
import logging
import string
import os.path
import os
import sys

from Ivy.version import __version__
import Ivy.utils

__program__ = 'ivy'
__author__ = 'Soh Ishiguro <yukke@g-language.org>'
__license__ = ''
__status__ = 'development'


class CommandLineParser(object):
    def __init__(self):
        desc = 'software package for identification of RNA editing sites based on massively parallel sequencing data'
        usage = 'usage: %prog [options]'
        fmt = IndentedHelpFormatter(indent_increment=2, max_help_position=60, width=120, short_first=1)
        self.parser = OptionParser(usage=usage, formatter=fmt, version=__version__, description=desc)
        #self.parse()
        
    def parse_basic_opt(self):
        # basic options
        self.parser.add_option('-f',
                               dest='fasta',
                               action='store',
                               metavar='',
                               nargs=1,
                               type='string',
                               help='Reference genome [fasta]'
                               )
        self.parser.add_option('-r',
                               dest='r_bams',
                               action='store',
                               metavar='',
                               type='string',
                               help='RNA-seq file(s) [bam]',
                               )
        self.parser.add_option('-b',
                               metavar='',
                               dest='d_bams',
                               action='store',
                               type='string',
                               help='DNA-seq file(s) [bam]',
                               )
        self.parser.add_option('-o',
                               dest='outname',
                               metavar='',
                               action='store',
                               type='string',
                               help='Output filename',
                               )
        self.parser.add_option('-l',
                               dest='regions',
                               action='store',
                               metavar='',
                               nargs=1,
                               type='string',
                               default='All',
                               help='Explore specify region [chr:start-end]',
                               # TODO:
                               # if not set -l, explore all region, default value is "all"?
                               )
        self.parser.add_option('-G',
                               metavar='',
                               dest='gtf',
                               action='store',
                               nargs=1,
                               type='string',
                               help='GTF file',
                               )
        self.parser.add_option('--one-based',
                               metavar='',
                               dest='one_based',
                               action='store_false',
                               default=False,
                               help='Genomic coordinate'
                               )
        self.parser.add_option('--num-threads',
                               metavar='',
                               dest='n_threads',
                               action='store',
                               nargs=1,
                               default=1,
                               type='int',
                               help='Number of threads [default: %default]',
                               )
        self.parser.add_option('--dry-run',
                               metavar='',
                               dest='dry_run',
                               action='store_true',
                               help='Dry run ivy'
                               )
        self.parser.add_option('--verbose',
                               metavar='',
                               action='store_true',
                               help='Show verbously messages'
                               )
        
    def parse_sample_opt(self):
        # sample options
        sample_group = OptionGroup(self.parser, 'Sample options')
        sample_group.add_option('--strand',
                                metavar='',
                                dest='strand',
                                action='store_true',
                                default=False,
                                help='Strand-specific seq. data is used. [default: %default]'
                                )
        sample_group.add_option('--ko-strain',
                                metavar='',
                                dest='ko_strain',
                                action='store_true',
                                default=False,
                                help='Adar null strain is used. [default: %default]'
                                )
        sample_group.add_option('--replicate',
                                metavar='',
                                dest='replicate',
                                action='store_false',
                                default=False,
                                help='Biological replicate is used [default: %default]'
                                )
        self.parser.add_option_group(sample_group)
        
    def parse_basic_filt_opt(self):
        # basic filter options
        basic_filter_group = OptionGroup(self.parser, 'Basic filter options')
        basic_filter_group.add_option('--min-ag-ratio',
                                      metavar='',
                                      dest='ag_ratio',
                                      action='store',
                                      nargs=1,
                                      default=0.1,
                                      type='float',
                                      help='Min A-to-G edit base ratio [default: %default]'
                                      )
        basic_filter_group.add_option('--min-rna-coverage',
                                      metavar='',
                                      dest='min_rna_cov',
                                      action='store',
                                      nargs=1,
                                      default=10,
                                      type='int',
                                      help='Min RNA read coverage [default: %default]'
                                      )
        basic_filter_group.add_option('--min-dna-coverage',
                                      metavar='',
                                      dest='min_dna_cov',
                                      action='store',
                                      nargs=1,
                                      default=20,
                                      type='int',
                                      help='Min DNA read coverage [default: %default]'
                                      )
        basic_filter_group.add_option('--rm-duplicated-read',
                                      metavar='',
                                      dest='is_duplicated',
                                      action='store_true',
                                      default=True,
                                      help='Remove duplicated reads [default: %default]'
                                      )
        basic_filter_group.add_option('--rm-deletion-read',
                                      metavar='',
                                      dest='is_deletion',
                                      action='store_true',
                                      default=True,
                                      help='Remove deletion reads [default: %default]'
                                      )
        basic_filter_group.add_option('--min-mapq',
                                      metavar='',
                                      dest='min_mapq',
                                      action='store',
                                      nargs=1,
                                      default=30,
                                      type='int',
                                      help='Min mapping quality [default: %default]'
                                      )
        basic_filter_group.add_option('--num-allow-type',
                                      metavar='',
                                      dest='num_type',
                                      action='store',
                                      nargs=1,
                                      default=1,
                                      type='int',
                                      help='Number of allowing base modification type [default: %default]'
                                      )
        basic_filter_group.add_option('--min-baq-rna',
                                      metavar='',
                                      dest='min_baq_r',
                                      action='store',
                                      nargs=1,
                                      default=28,
                                      type='int',
                                      help='Min base call quality in RNA [default: %default]'
                                      )
        basic_filter_group.add_option('--min-baq-dna',
                                      metavar='',
                                      dest='min_baq_d',
                                      action='store',
                                      nargs=1,
                                      default=28,
                                      type='int',
                                      help='Min base call quality in DNA [default: %default]'
                                      )
        self.parser.add_option_group(basic_filter_group)
        
    def parse_stat_filt_opt(self):
        # statistical filters options
        stat_filter_group = OptionGroup(self.parser, 'Statistical filter options')
        stat_filter_group.add_option('--sig-level',
                                     metavar='',
                                     dest='sig_level',
                                     action='store',
                                     default=0.05,
                                     nargs=1,
                                     type='float',
                                     help='Significance level [default: %default]'
                                     )
        stat_filter_group.add_option('--base-call-bias',
                                     metavar='',
                                     dest='baq_bias',
                                     action='store_true',
                                     default=False,
                                     help='Consider base call bias [default: %default]'
                                     )
        stat_filter_group.add_option('--strand-bias',
                                     metavar='',
                                     dest='strand_bias',
                                     action='store_true',
                                     default=False,
                                     help='Consider strand bias [default: %default]'
                                     )
        stat_filter_group.add_option('--positional-bias',
                                     metavar='',
                                     dest='pos_bias',
                                     action='store_true',
                                     default=False,
                                     help='Consider positional bias [default: %default]'
                                     )
        self.parser.add_option_group(stat_filter_group)
        
    def parse_ext_filt_opt(self):
        # extended options
        ext_filter_group = OptionGroup(self.parser, 'Extended filter options')
        ext_filter_group.add_option('--blat-collection',
                                    metavar='',
                                    dest='blat',
                                    action='store_false',
                                    default=False,
                                    help='Reduce mis-alignment with Blat [default: %default]'
                                    )
        ext_filter_group.add_option('--snp',
                                    metavar='',
                                    dest='snp_file',
                                    action='store',
                                    help='Exclude variation sites [vcf]'
                                    )
        ext_filter_group.add_option('--ss-num',
                                    metavar='',
                                    dest='ss_num',
                                    action='store',
                                    nargs=1,
                                    default=5,
                                    type='int',
                                    help='Exclude site around the splice sistes [default: %defaultbp]'
                                    )
        ext_filter_group.add_option('--trim-n',
                                    metavar='',
                                    dest='trim_n',
                                    action='store',
                                    nargs=1,
                                    default=10,
                                    type='int',
                                    help='Do not call Nbp in up/down read [default: %defaultbp]'
                                    )
        ext_filter_group.add_option('--mask-repeat',
                                    metavar='',
                                    dest='is_mask_rep',
                                    action='store_false',
                                    default=False,
                                    help='Mask repeat sequence [default: %default]'
                                    )
        self.parser.add_option_group(ext_filter_group)

    def ivy_parse_options(self):
        self.parse_basic_opt()
        self.parse_ext_filt_opt()
        self.parse_sample_opt()
        self.parse_basic_filt_opt()
        self.parse_stat_filt_opt()
        
        (opt, args) = self.parser.parse_args()

        # enable --dry-run or not
        if opt.dry_run:
            print "### All options with values ###"
            for k, v in self.parser.values.__dict__.iteritems():
                print '[' + k + ']:', v
            die()

        #############################
        ### Check required params ###
        #############################
        passed_params = {}
        
        # fasta file, -f
        if not opt.fasta:
            self.parser.error('[-f] Reference fasta file is a required argument')
        elif self._ok_file(opt.fasta):
            passed_params.update({'fasta': opt.fasta})
        elif not self._ok_file(opt.fasta):
            self.parser.error(opt.fasta + " is not found or writable file!")
            
        # RNA-seq bam file, r_bams
        if not opt.r_bams:
            self.parser.error('[-r] RNA-seq bam file is a required argument')
        elif self._ok_file(opt.r_bams):
            passed_params.update({'r_bams': opt.r_bams})
        elif not self._ok_file(opt.r_bams):
            self.parser.error(opt.r_bams + " is not found or writable file!")
        
        # output filename, -o
        if opt.outname:
            passed_params.update({'outname': opt.outname})
            #self.parser.error('[-o] Output filename is a required argument')
        else:
            default_filename = 'ivy_run.log'
            passed_params.update({'outname': default_filename})
        
        ###########################
        ### Check basic options ###
        ###########################
        # -l, regions
        if opt.regions == 'All':
            # default value: all
            passed_params.update({'region': opt.regions})
        elif self._is_region(opt.regions) is not None:
            # specified region: e.g. chr1:1-1000
            passed_params.update({'region': self._is_region(opt.regions)})
        else:
            # error
            self.parser.error("faild to set region")

        # gtf file, -G
        if opt.gtf:
            passed_params.update({'gtf': opt.gtf})
        else:
            passed_params.update({'gtf': None})

        # --one_based
        if opt.one_based is True:
            passed_params.update({'one_based': opt.one_based})
        elif opt.one_based is False:
            passed_params.update({'one_based': False})

        # --num-threads
        if opt.n_threads:
            passed_params.update({'n_threads': opt.n_threads})

        ############################
        ### Check sample options ###
        ###########################
        # --strand
        if opt.strand is True:
            passed_params.update({'strand': opt.strand})
        elif opt.strand is False:
            passed_params.update({'strand': False})

        # --ko-strain
        if opt.ko_strain is True:
            passed_params.update({'strand': opt.strand})
        elif opt.ko_strain is False:
            passed_params.update({'strand': False})

        # --replicate
        if opt.replicate is True:
            passed_params.update({'replicate': opt.replicate})
        elif opt.replicate is False:
            passed_params.update({'replicate': opt.replicate})

        ############################
        ### Basic filter options ###
        ############################
        # --min-ag-ratio
        if opt.ag_ratio:
            passed_params.update({'ag_ratio': opt.ag_ratio})

        # --min-rna-coverage
        if opt.min_rna_cov:
            passed_params.update({'min_rna_cov': opt.min_rna_cov})

        # --min_dna_coverage
        if opt.min_dna_cov:
            passed_params.update({'min_dna_cov': opt.min_dna_cov})

        # --rm-duplicated-read
        if opt.is_duplicated:
            passed_params.update({'is_duplicated': opt.is_duplicated})

        # --rm-deletion-read
        if opt.is_deletion:
            passed_params.update({'is_deletion': opt.is_deletion})

        # --min-mapq
        if opt.min_mapq:
            passed_params.update({'min_mapq': opt.min_mapq})

        # --num-allow-type
        if opt.num_type:
            passed_params.update({'num_type': opt.num_type})

        # --min-baq-rna
        if opt.min_baq_r:
            passed_params.update({'min_baq_rna': opt.min_baq_r})

        # --min-baq-dna
        if opt.min_baq_d:
            passed_params.update({'min_baq_dba': opt.min_baq_d})

        ##################################
        ### Statistical filter options ###
        ##################################
        # --sig-level
        if opt.sig_level:
            passed_params.update({'sig_level': opt.sig_level})
            
        # base-call-bias
        if opt.baq_bias:
            passed_params.update({'baq_bias': opt.baq_bias})

        # strand-bias
        if opt.strand_bias:
            passed_params.update({'strnad_bias': opt.strand_bias})
        
        # --potitional-bias
        if opt.pos_bias:
            passed_params.update({'has bias': opt.pos_bias})

        ###########################
        ### Ext. filter options ###
        ###########################
        # --blat-collection
        if opt.blat:
            passed_params.update({'blat': opt.blat})
        elif opt.blat is False:
            passed_params.update({'blat': False})

        # --snp
        if opt.snp_file:
            passed_params.update({'snp_file': opt.snp})
            
        # --ss-num
        if opt.ss_num:
            passed_params.update({'ss_num': opt.ss_num})

        # --trim-n
        if opt.trim_n:
            passed_params.update({'trim_n': opt.trim_n})

        # --mask-repeat
        if opt.is_mask_rep:
            passed_params.update({'is_mask_rep': opt.is_mask_rep})
 
        return passed_params

    def _ok_file(self, filename):
        if os.path.isfile(filename) and os.access(filename, os.R_OK):
            return True
        else:
            return False

    def _is_region(self, regions):
        # TODO: needed to unittest!
        if regions == 'All':
            return True
        elif regions:
            try:
                (chrom, pos) = regions.split(':')
            except ValueError:
                self.parser.error("[" + regions+ "]" + ' is lacked chromosome or position value')
            
            if not chrom:
                self.parser.error(regions + 'is invalid chromosome name')
                return False
            elif not pos:
                self.parser.error(regions + 'is invalid position')
                return False
            else:
                try:
                    (start, end) = pos.split('-')
                except ValueError:
                    self.parser.error("position must be \'start-end\'")
                
                if start.isdigit() and end.isdigit():
                    if start < end:
                        # everything is fine
                        return {'chrom': str(chrom), 'start': int(start), 'end': int(end)}
                    elif start > end:
                        self.parser.error('end:' + end + ' is greater than ' + 'start:' + start)
                        return False
                    elif start == end:
                        self.parser.error("start:" + start + ", end:" + end + " is same values")
                        return False
                else:
                    self.parser.error(regions + ' in pos is not numetric (expected integer)')
                    return False
        else:
            return False
            
def die(msg=''):
    raise SystemExit(msg)
