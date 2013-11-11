from Ivy.version import __version__
import Ivy.utils
from optparse import OptionParser, OptionGroup, HelpFormatter, IndentedHelpFormatter
import os.path

__program__ = 'ivy'
__author__ = 'Soh Ishiguro <yukke@g-language.org>'
__license__ = ''
__status__ = 'development'


class MyParser(OptionParser):
    def format_epilog(self, formatter):
        return self.epilog
        
def run():
    desc = "software package for identification of RNA editing sites based on massively parallel sequencing data\n"
    usage = '%prog'

    fmt = IndentedHelpFormatter(indent_increment=2, max_help_position=60, width=120, short_first=1)
    parser = OptionParser(usage=usage, formatter=fmt)
    basic_filter_group = OptionGroup(parser, 'Basic filter options')
    ext_filter_group = OptionGroup(parser, 'Extended filter options')
    stat_filter_group = OptionGroup(parser, 'Statistical filter options')
    sample_group = OptionGroup(parser, 'Sample options')

    # Basic options
    parser.add_option('-f',
                      dest='fasta',
                      action='store',
                      metavar='',
                      nargs=1,
                      help='set reference genome [fasta]'
                      )
    parser.add_option('-r',
                      dest='r_bams',
                      action='store',
                      metavar='',
                      nargs='+',
                      help='RNA-seq file(s) [bam]',
                      )
    parser.add_option('-b',
                      metavar='',
                      dest='d_bams',
                      action='store',
                      nargs='+',
                      help='DNA-seq file(s) [bam]',
                      )
    parser.add_option('-l',
                      dest='regions',
                      action='store',
                      metavar='',
                      nargs=2,
                      help='Explore specify region [chr:start chr:end]'
                      )
    parser.add_option('-G',
                      metavar='',
                      dest='gtf',
                      action='store',
                      nargs=1,
                      help='GTF file',
                      )
    parser.add_option('--num_threads',
                      metavar='',
                      dest='n_threads',
                      action='store',
                      nargs=1,
                      default=1,
                      help='Number of threads [default: %default]',
                      )
    parser.add_option('--verbose',
                      metavar='',
                      dest='verbose',
                      default=False,
                      help='Show verbously messages'
                      )
    
    # sample options
    sample_group.add_option('--strand',
                            metavar='',
                            action='store',
                            nargs=1,
                            default=False,
                            help='Strand-specific seq. data is used. [default: %default]'
                            )
    sample_group.add_option('--ko_strain',
                            metavar='',
                            action='store',
                            nargs=1,
                            default=False,
                            help='Adar null strain is used. [default: %default]'
                            )
    sample_group.add_option('--replicate',
                            metavar='',
                            action='store',
                            nargs=1,
                            default=False,
                            help='Biological replicate is used [default: %default]'
                            )

    # basic filter options
    basic_filter_group.add_option('--min_ag_ratio',
                                  metavar='',
                                  dest='ag_ratio',
                                  action='store',
                                  nargs=1,
                                  default=0.1,
                                  help='Min A-to-G edit base ratio [default: %default]'
                                  )
    basic_filter_group.add_option('--min_rna_coverage',
                                  metavar='',
                                  dest='min_rna_cov',
                                  action='store',
                                  nargs=1,
                                  default=10,
                                  help='Min RNA read coverage [default: %default]'
                                  )
    basic_filter_group.add_option('--min_dna_coverage',
                                  metavar='',
                                  dest='min_dna_cov',
                                  action='store',
                                  nargs=1,
                                  default=20,
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
                                  help='Min mapping quality [default: %default]'
                                  )
    basic_filter_group.add_option('--num_allow_type',
                                  metavar='',
                                  dest='num_type',
                                  action='store',
                                  nargs=1,
                                  default=1,
                                  help='Number of allowing base modification type [default: %default]'
                                  )
    basic_filter_group.add_option('--min_baq_rna',
                                  metavar='',
                                  dest='min_baq_r',
                                  action='store',
                                  nargs=1,
                                  default=28,
                                  help='Min base call quality in RNA [default: %default]'
                                  )
    basic_filter_group.add_option('--min_baq_dna',
                                  metavar='',
                                  dest='min_baq_d',
                                  action='store',
                                  nargs=1,
                                  default=28,
                                  help='Min base call quality in DNA [default: %default]'
                                  )
    stat_filter_group.add_option('--sig_level',
                                 metavar='',
                                 dest='sig_level',
                                 action='store',
                                 nargs=1,
                                 default=0.05,
                                 help='Significance level [default: %default]'
                                 )
    stat_filter_group.add_option('--base_call_bias',
                                 metavar='',
                                 dest='baq_bias',
                                 action='store',
                                 nargs=1,
                                 default=True,
                                 help='Consider base call bias [default: %default]'
                                 )
    stat_filter_group.add_option('--strand_bias',
                                 metavar='',
                                 dest='strand_bias',
                                 action='store',
                                 nargs=1,
                                 default=True,
                                 help='Consider strand bias [default: %default]'
                                 )
    stat_filter_group.add_option('--positional_bias',
                                 metavar='',
                                 dest='pos_bias',
                                 action='store',
                                 nargs=1,
                                 default=True,
                                 help='Consider positional bias [default: %default]'
                                 )
    ext_filter_group.add_option('--blat_collection',
                                metavar='',
                                dest='blat',
                                action='store',
                                nargs=1,
                                default=False,
                                help='Reduce mis-alignment with Blat [default: %default]'
                                )
    ext_filter_group.add_option('--snp',
                                metavar='',
                                dest='snp',
                                action='store',
                                nargs=1,
                                default=False,
                                help='Exclude sites within SNP locations [default: %default]'
                                )
    ext_filter_group.add_option('--ss_num',
                                metavar='',
                                dest='ss_num',
                                action='store',
                                nargs=1,
                                default=5,
                                help='Exclude site around the splice sistes [default: %defaultbp]'
                                )
    ext_filter_group.add_option('--trim_n',
                                metavar='',
                                dest='trim_n',
                                action='store',
                                nargs=1,
                                default=10,
                                help='Do not call Nbp in up/down read [default: %defaultbp]'
                                )
    ext_filter_group.add_option('--mask_repeat',
                                metavar='',
                                dest='is_mask_rep',
                                action='store',
                                nargs=1,
                                default=False,
                                help='Mask repeat sequence [default: %default]'
                                )

    
    parser.add_option_group(basic_filter_group)
    parser.add_option_group(stat_filter_group)
    parser.add_option_group(ext_filter_group)
    parser.add_option_group(sample_group)
    (options, args) = parser.parse_args()    

if __name__ == '__main__':
    run()
