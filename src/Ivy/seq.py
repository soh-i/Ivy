import string
import re
import os.path

class Fasta(object):
    def __init__(self, filename):
        if os.path.isfile(filename):
            self.filename = filename
        else:
            raise RuntimeError("{fa:s} isnot found".format(fa=fastafile))

    
        
