from Ivy.version import __version__
from Ivy.alignment.stream import RNASeqAlignmentStream, DNASeqAlignmentStream
from Ivy.commandline.parse_ivy_opts import CommandLineParser
from Ivy.annotation.writer import VCFWriteHeader
from pprint import pprint as p
import logging

__program__ = 'run_ivy'
__author__ = 'Soh Ishiguro <yukke@g-language.org>'
__license__ = ''
__status__ = 'development'
__version__ = __version__

def run():
    logger = logging.getLogger(__name__)    
    parse = CommandLineParser()
    params = parse.ivy_parse_options()
    vcf = VCFWriteHeader(params)
    #vcf.make_vcf_header()

    logger.debug("Beginning Ivy run (v." + __version__ + ")" )
    if params.r_bams:
        logger.debug("Loading RNA-seq bam file '{0}'".format(params.r_bams))
        rna_pileup_alignment = RNASeqAlignmentStream(params)
        for rna in rna_pileup_alignment.filter_stream():
            print to_tab(rna)
            #print to_tsv(rna)

    if params.d_bams:
        logger.debug("Loading DNA-seq bam file '{0}'".format(params.d_bams))
        dna_pileup_alignment = DNASeqAlignmentStream(params)
        for dna in dna_pileup_alignment.filter_stream():
            #pprint(dna)
            print to_tsv(dna)

    #with open(params.outname, 'w') as f:
    #f.write()

    
def pprint(data, *args, **kwargs):
    '''
    Simple pretty print yielded pileuped data
    '''
    
    print '{'
    for key in data:
        print '  {key}({types}): {val}'.format(
            key=key, val=data[key], types=type(data[key]).__name__)
    print '}'

def to_tsv(data, *args, **kwargs):
    entory = ''
    for key in data:
        if isinstance(data[key], str):
            entory += data[key]+"\t"
        elif isinstance(data[key], int) or isinstance(data[key], float):
            entory += str(data[key])+"\t"
        elif isinstance(data[key], tuple):
            entory += ','.join([str(_) for _ in data[key]]) + "\t"
    return entory
    
            
def to_tab(data, *args, **kwargs):
    return '{chrom:}\t{pos:}\t{ref:}\t{alt:}\t{dp4:}'.format(
        chrom=data.get('chrom'),
        pos=data.get('pos'),
        ref=data.get('ref'),
        alt=data.get('alt'),
        dp4=",".join([str(_) for _ in data.get('dp4')]))
    
            
