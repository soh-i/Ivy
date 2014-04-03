from argparse import ArgumentParser
import os.path
import sys
from Ivy.version import __version__
from Ivy.analysis_settings import EDIT_BENCH_SETTINGS

__program__ = 'ivy_benchmark'
__author__ = 'Soh Ishiguro <yukke@g-language.org>'
__license__ = ''
__status__ = 'development'

# Do not call SystemExit, instead of throwing exception
class ArgumentParserError(Exception): pass
    
class ThrowingArgumentParser(ArgumentParser):
    def error(self, message):
        raise ArgumentParserError(message)

        
def parse_bench_opts():
    desc = "Benchmarking test for detected RNA editing sites based on HTSeq data to evaluate detection params."
    parser = ThrowingArgumentParser(description=desc,
                            prog=__program__,
                            )
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--vcf',
                       dest='vcf_file',
                       action='store',
                       nargs='+',
                       metavar='',
                       help='VCF file(s).')
    group.add_argument('--csv',
                       dest='csv_file',
                       action='store',
                       nargs='+',
                       metavar='',
                       help='CSV file(s), for ***debug mode***.')
    parser.add_argument('--source',
                        required=False,
                        default=EDIT_BENCH_SETTINGS['APP']['SOURCE'],
                        dest='source',
                        action='store',
                        metavar='',
                        help='To use specific sample/tissue/cell line. [default: {0}]'.format(
                            EDIT_BENCH_SETTINGS['APP']['SOURCE']))
    parser.add_argument('--sp',
                        required=True,
                        default=EDIT_BENCH_SETTINGS['APP']['SP'],
                        dest='sp',
                        action='store',
                        metavar='species',
                        help='Species + genome version. (eg. human_hg19)')
    parser.add_argument('--plot',
                        required=False,
                        default=EDIT_BENCH_SETTINGS['APP']['PLOT'],
                        action='store_true',
                        help='Make a precision-recall plot. [default: {0}]'.format(
                            EDIT_BENCH_SETTINGS['APP']['PLOT']))
    parser.add_argument('--out',
                        dest='out',
                        default=EDIT_BENCH_SETTINGS['APP']['OUT'],
                        required=False,
                        action='store',
                        metavar='out',
                        help='Output file name. [default: {0}]'.format(
                            EDIT_BENCH_SETTINGS['APP']['OUT']))
    parser.add_argument('--version',
                        action='version',
                        help='Show program version number and exit.',
                        version=__version__)
    return parser.parse_args()
