from Ivy.version import __version__
import Ivy.utils
from optparse import OptionParser, OptionGroup
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
    
    parser = OptionParser(usage)
    input_group = OptionGroup(parser, 'Options')
    basic_filter_group = OptionGroup(parser, 'Basic filter options')
    stat_filter_group = OptionGroup(parser, 'Statistical filter options')

    # Input params
    input_group.add_option('-f', '--fasta',
                           dest='fasta',
                           action='store',
                           metavar=' ',
                           nargs=1,
                           help='set reference genome [fasta]'
                           )
    input_group.add_option('-r', '--rbam',
                           dest='r_bams',
                           action='store',
                           metavar=' ',
                           nargs='+',
                           help='RNA-seq file(s) [bam]',
                           )
    input_group.add_option('-d', '--dbam',
                           metavar=' ',
                           dest='d_bams',
                           action='store',
                           nargs='+',
                           help='DNA-seq file(s) [bam]',
                           )
    input_group.add_option('-R', '--region',
                           dest='regions',
                           action='store',
                           metavar=' ',
                           nargs=2,
                           help='Explore specify region [chr:start chr:end]'
                           )
    input_group.add_option('-G', '--GTF',
                           metavar=' ',
                           dest='gtf',
                           action='store',
                           nargs=1,
                           help='GTF',
                           )
    input_group.add_option('-p', '--num_threads',
                           metavar=' ',
                           dest='n_threads',
                           action='store',
                           nargs=1,
                           help='Number of threads',
                           )
    basic_filter_group.add_option('-e', '--max_edit_ratio',
                                  metavar=' ',
                                  dest='edit_ratio',
                                  action='store',
                                  nargs=1,
                                  help='Max edit base ratio'
                                  )
    
    parser.add_option_group(input_group)
    parser.add_option_group(basic_filter_group)
    (options, args) = parser.parse_args()    

if __name__ == '__main__':
    run()
