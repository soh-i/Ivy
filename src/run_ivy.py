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

    # Optiopns
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
                      help='GTF',
                      )
    parser.add_option('--num_threads',
                      metavar='',
                      dest='n_threads',
                      action='store',
                      nargs=1,
                      help='Number of threads',
                      )
    parser.add_option('--verbose',
                      metavar='',
                      dest='verbose',
                      default=False,
                      help='Show verbously messages'
                      )
                      
    basic_filter_group.add_option('--max_edit_ratio',
                                  metavar='',
                                  dest='edit_ratio',
                                  action='store',
                                  nargs=1,
                                  help='Max edit base ratio'
                                  )
    basic_filter_group.add_option('--min_coverage',
                                  metavar='',
                                  dest='min_coverage',
                                  action='store',
                                  nargs=1,
                                  help='Min read coverage'
                                  )
    stat_filter_group.add_option('--base_call_bias',
                                 metavar='',
                                 dest='baq_bias',
                                 action='store',
                                 nargs=1,
                                 help='Consider base call bias'
                                 )
    stat_filter_group.add_option('--strand_bias',
                                 metavar='',
                                 dest='strand_bias',
                                 action='store',
                                 nargs=1,
                                 help='Consider strand bias'
                                 )
    ext_filter_group.add_option('--blat_collection',
                                metavar='',
                                dest='blat',
                                action='store',
                                nargs=1,
                                help='Reduce mis-alignment with blat'
                                )
    ext_filter_group.add_option('--snp',
                                metavar='',
                                dest='snp',
                                action='store',
                                nargs=1,
                                help='Exclude SNP locations'
                                )
    
    parser.add_option_group(basic_filter_group)
    parser.add_option_group(stat_filter_group)
    parser.add_option_group(ext_filter_group)
    (options, args) = parser.parse_args()    

if __name__ == '__main__':
    run()
