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
    parser = argparse.ArgumentParser(description=desc,
                                     prog=__program__,
                                     )
    parser.add_argument('--version', action='version', version=__version__)
    parser.add_argument('--ref', required=True, metavar='', dest='ref', action='store', help='set reference genome [fasta]')
    parser.add_argument('--bam', required=True, metavar='', dest='bams', action='store', nargs='+', help='set BAM file(s)')
    
    args = parser.parse_args()


if __name__ == '__main__':
    run()
