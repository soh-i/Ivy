from Ivy.version import __version__
from Ivy.alignment.stream import AlignmentConfig
from Ivy.parse_opt import CommandLineParser
import logging

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
        self.parse()

        #align_conf = AlignmentConfig(passed_params)
        
    
