from Ivy.version import __version__
from Ivy.alignment.stream import (
    AlignmentConfig,
    AlignmentStream,
    AlignmentPreparation,
    )
from Ivy.parse_opt import CommandLineParser
import logging
import pprint

__program__ = 'run_ivy'
__author__ = 'Soh Ishiguro <yukke@g-language.org>'
__license__ = ''
__status__ = 'development'


class Ivy(CommandLineParser):
    def __init__(self):
        logging.basicConfig(level=logging.ERROR, format="%(asctime)s %(message)s")
        logging.error('Job started')

    def run(self):
        CommandLineParser.__init__(self)
        params = self.ivy_parse_options()

        # just test code is here
        params.update({'chrom': '21', 'start': 47721030, 'end': 47721057, 'one_based': True})
        align_conf = AlignmentConfig(params)
        print align_conf.cl_params
        
        stream = AlignmentStream(align_conf.cl_params)
        pp = pprint.PrettyPrinter(indent=6)
        for rna in stream.pileup_stream():
            pp.pprint(rna)
         

        
    
