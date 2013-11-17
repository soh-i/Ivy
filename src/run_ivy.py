from Ivy.version import __version__
from Ivy.alignment.stream import AlignmentStream
from Ivy.parse_opt import CommandLineParser
import logging
import pprint

__program__ = 'run_ivy'
__author__ = 'Soh Ishiguro <yukke@g-language.org>'
__license__ = ''
__status__ = 'development'


class Ivy(object):
    def __init__(self):
        logging.basicConfig(level=logging.ERROR, format="%(asctime)s %(message)s")
        logging.error('[Ivy] Job started')
        
    def run(self):
        pp = pprint.PrettyPrinter(indent=6)
        parse = CommandLineParser()
        params = parse.ivy_parse_options()
        stream = AlignmentStream(params)
                
        for rna in stream.pileup_stream():
            print "%s\t%d\t%s\t%s\t%d\t%f" % (rna['chrom'], rna["pos"], rna['ref'], rna['alt'], rna['coverage'], rna['mismatch_ratio'])
