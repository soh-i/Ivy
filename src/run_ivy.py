from Ivy.version import __version__
import Ivy.utils
import argparse
import os.path

__program__ = 'ivy'
__author__ = 'Soh Ishiguro <yukke@g-language.org>'
__license__ = ''
__status__ = 'development'


def run():
    desc = "software package for identification of RNA editing sites based on massively parallel sequencing data."
    parser = argparse.ArgumentParser(description=desc, prog=__program__)

    # Basic param
    parser.add_argument('--version', action='version', version=__version__)
    parser.add_argument('-f', '--fasta', required=True, metavar='', dest='ref', action='store', nargs=1, help='set reference genome [fasta]')
    parser.add_argument('-r', '--rbam', required=True, metavar='', dest='r_bams', action='store', nargs='+', help='RNA-seq file(s) [bam]')
    parser.add_argument('-d', '--dbam', required=True, metavar='', dest='d_bams', action='store', nargs='+', help='DNA-seq file(s) [bam]')
    parser.add_argument('-l', '--region', required=False, metavar='', dest='region', action='store', nargs=1, help='Explore specify region')
    parser.add_argument('-G', '--GTF', required=False, metavar='', dest='gtf', action='store', nargs=1, help='GTF for annotation')
    parser.add_argument('-p', '--num_threads', required=False, metavar='', dest='n_threads', action='store', nargs=1, help='Number of threads')
    parser.add_argument('-o', '--out', required=True, metavar='', dest='outname', action='store', nargs=1, help='Output filename')
    
    parser.add_argument('-o', '--out', required=True, metavar='', dest='outname', action='store', nargs=1, help='Output filename')
    parser.add_argument('-o', '--out', required=True, metavar='', dest='outname', action='store', nargs=1, help='Output filename')
    parser.add_argument('-o', '--out', required=True, metavar='', dest='outname', action='store', nargs=1, help='Output filename')
    
    
    
    args = parser.parse_args()


if __name__ == '__main__':
    run()
