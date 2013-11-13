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
                               help='Explore specify region [chr:start-end]'
                               )
        self.parser.add_option('-G',
                               metavar='',
                               dest='gtf',
                               action='store',
                               nargs=1,
                               type='string',
                               help='GTF file',
                               )
        self.parser.add_option('--one_based',
                               metavar='',
                               action='store_false',
                               help='Genomic coordinate'
                               )
        self.parser.add_option('--num_threads',
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
                               action='store_false',
                               help='Show verbously messages'
                               )
        
    def parse_sample_opt(self):
        # sample options
        sample_group = OptionGroup(self.parser, 'Sample options')
        sample_group.add_option('--strand',
                                metavar='',
                                default=False,
                                action='store_false',
                                help='Strand-specific seq. data is used. [default: %default]'
                                )
        sample_group.add_option('--ko_strain',
                                metavar='',
                                default=False,
                                action='store_false',
                                help='Adar null strain is used. [default: %default]'
                                )
        sample_group.add_option('--replicate',
                                metavar='',
                                default=False,
                                action='store_false',
                                help='Biological replicate is used [default: %default]'
                                )
        self.parser.add_option_group(sample_group)
        
    def parse_basic_filt_opt(self):
        # basic filter options
        basic_filter_group = OptionGroup(self.parser, 'Basic filter options')
        basic_filter_group.add_option('--min_ag_ratio',
                                      metavar='',
                                      dest='ag_ratio',
                                      action='store',
                                      nargs=1,
                                      default=0.1,
                                      type='float',
                                      help='Min A-to-G edit base ratio [default: %default]'
                                      )
        basic_filter_group.add_option('--min_rna_coverage',
                                      metavar='',
                                      dest='min_rna_cov',
                                      action='store',
                                      nargs=1,
                                      default=10,
                                      type='int',
                                      help='Min RNA read coverage [default: %default]'
                                      )
        basic_filter_group.add_option('--min_dna_coverage',
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
                                      action='store',
                                      nargs=1,
                                      default=True,
                                      help='Remove duplicated reads [default: %default]'
                                      )
        basic_filter_group.add_option('--rm-deletion-read',
                                      metavar='',
                                      dest='is_deletion',
                                      action='store',
                                      nargs=1,
                                      default=True,
                                      help='Remove deletion reads [default: %default]'
                                      )
        basic_filter_group.add_option('--min_mapq',
                                      metavar='',
                                      dest='min_mapq',
                                      action='store',
                                      nargs=1,
                                      default=30,
                                      type='int',
                                      help='Min mapping quality [default: %default]'
                                      )
        basic_filter_group.add_option('--num_allow_type',
                                      metavar='',
                                      dest='num_type',
                                      action='store',
                                      nargs=1,
                                      default=1,
                                      type='int',
                                      help='Number of allowing base modification type [default: %default]'
                                      )
        basic_filter_group.add_option('--min_baq_rna',
                                      metavar='',
                                      dest='min_baq_r',
                                      action='store',
                                      nargs=1,
                                      default=28,
                                      type='int',
                                      help='Min base call quality in RNA [default: %default]'
                                      )
        basic_filter_group.add_option('--min_baq_dna',
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
        stat_filter_group.add_option('--sig_level',
                                     metavar='',
                                     dest='sig_level',
                                     action='store',
                                     default=0.05,
                                     nargs=1,
                                     type='float',
                                     help='Significance level [default: %default]'
                                     )
        stat_filter_group.add_option('--base_call_bias',
                                     metavar='',
                                     dest='baq_bias',
                                     action='store_true',
                                     default=True,
                                     help='Consider base call bias [default: %default]'
                                     )
        stat_filter_group.add_option('--strand_bias',
                                     metavar='',
                                     dest='strand_bias',
                                     action='store_true',
                                     default=True,
                                     help='Consider strand bias [default: %default]'
                                     )
        stat_filter_group.add_option('--positional_bias',
                                     metavar='',
                                     dest='pos_bias',
                                     action='store_true',
                                     default=True,
                                     help='Consider positional bias [default: %default]'
                                     )
        self.parser.add_option_group(stat_filter_group)
        
    def parse_ext_filt_opt(self):
        # extended options
        ext_filter_group = OptionGroup(self.parser, 'Extended filter options')
        ext_filter_group.add_option('--blat_collection',
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
        ext_filter_group.add_option('--ss_num',
                                    metavar='',
                                    dest='ss_num',
                                    action='store',
                                    nargs=1,
                                    default=5,
                                    type='int',
                                    help='Exclude site around the splice sistes [default: %defaultbp]'
                                    )
        ext_filter_group.add_option('--trim_n',
                                    metavar='',
                                    dest='trim_n',
                                    action='store',
                                    nargs=1,
                                    default=10,
                                    type='int',
                                    help='Do not call Nbp in up/down read [default: %defaultbp]'
                                    )
        ext_filter_group.add_option('--mask_repeat',
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
                print k+':', v
            die()

        ###
        ### Check input some files ###
        ###
        passed_params = {}
        
        # fasta file
        if not opt.fasta:
            self.parser.error('[-f] Reference fasta file is a required argument')
        elif self._ok_file(opt.fasta):
            passed_params.update({'fasta': opt.fasta})
        elif not self._ok_file(opt.fasta):
            self.parser.error(opt.fasta + " is not found or writable file!")
            
        # bam file
        if not opt.r_bams:
            self.parser.error('[-r] RNA-seq bam file is a required argument')
        elif self._ok_file(opt.r_bams):
            passed_params.update({'r_bams': opt.r_bams})
        elif not self._ok_file(opt.r_bams):
            self.parser.error(opt.r_bams + " is not found or writable file!")
        
        # output filename
        if opt.outname:
            passed_params.update({'outname': opt.outname})
            #self.parser.error('[-o] Output filename is a required argument')
        else:
            default_filename = 'ivy_run.log'
            passed_params.update({'outname': default_filename})

        ###
        ### Check basic opts
        ###

        # -l, regions
        if opt.regions:
            if self._is_region(opt.regions):
                passed_params.update({'region':opt.regions})
                print passed_params
                die("OK")
        
        return passed_params

    def _ok_file(self, filename):
        if os.path.isfile(filename) and os.access(filename, os.R_OK):
            return True
        else:
            return False

    def _is_region(self, regions):
        # TODO: needed to unittest!
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
                    return True
                elif start > end:
                    self.parser.error('end:' + end + ' is greater than ' + 'start:' + start)
                    return False
                elif start == end:
                    self.parser.error("start:" + start + ", end:" + end + " is same values")
                    return False
            else:
                self.parser.error(regions + ' in pos is not numetric (expected integer)')
                return False
            
def die(msg=''):
    raise SystemExit(msg)
