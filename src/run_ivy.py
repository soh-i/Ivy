from Ivy.version import __version__
from Ivy.alignment.stream import AlignmentStream
from Ivy.alignment.filter import AlignmentReadFilter
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
        params = CommandLineParser()
        stream = AlignmentStream(params.ivy_parse_options())
        #stream = AlignmentReadFilter(params.ivy_parse_options())
        
        for rna in stream.pileup_stream():
            #pp.pprint(rna)
            if rna['mismatches'] > 10:
                print "%s\t%s\t%s\t%s" % (rna['chrom'], rna["pos"], rna['ref'], rna['alt'])
                #raise SystemExit()
            
