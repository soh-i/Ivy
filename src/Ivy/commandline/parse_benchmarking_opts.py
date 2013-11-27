from argparse import ArgumentParser
import os.path
import sys
from Ivy.version import __version__

__program__ = 'ivy_benchmark'
__author__ = 'Soh Ishiguro <yukke@g-language.org>'
__license__ = ''
__status__ = 'development'


def parse_bench_opts():
    '''
    Parse command line arguments from ivy_bench.
    
    Returns:
     argparse object
    
    Examples:
     >>> args = parse_bench_opts()
    '''
    
    desc = "Benchmarking test for detected RNA editing sites based on HTSeq data to evaluate detection params."
    
    parser = ArgumentParser(description=desc,
                            prog=__program__,
                            )
    group = parser.add_mutually_exclusive_group(required=True)
    
    group.add_argument('--vcf',
                       dest='vcf_file',
                       action='store',
                       nargs='+',
                       metavar='',
                       help='VCF file(s)',
                       )
    group.add_argument('--csv',
                       dest='csv_file',
                       action='store',
                       nargs='+',
                       metavar='',
                       help='CSV file(s), For debug mode',
                       )
    parser.add_argument('--source',
                        required=False,
                        default='All',
                        dest='source',
                        action='store',
                        metavar='',
                        help='specific sample/tissue/cell line [default: All]',
                        )
    parser.add_argument('--sp',
                        required=True,
                        dest='sp',
                        metavar='species',
                        action='store',
                        help='species and genome version (eg. human_hg19)',
                        )
    parser.add_argument('--plot',
                        required=False,
                        default=False,
                        action='store_true',
                        help='plot benchmarking stats [default: Off]',
                        )
    parser.add_argument('--out',
                        dest='out',
                        required=False,
                        action='store',
                        metavar='out',
                        help='output name',
                        )
    parser.add_argument('--version',
                        action='version',
                        version=__version__
                        )
    return parser.parse_args()
    
