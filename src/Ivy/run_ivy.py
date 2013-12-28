from Ivy.alignment.stream import RNASeqAlignmentStream, DNASeqAlignmentStream
from Ivy.commandline.parse_ivy_opts import CommandLineParser
from Ivy.annotation.writer import VCFWriteHeader
from Ivy.seq import Fasta
from Ivy.version import __version__
from pprint import pprint as p
from multiprocessing import Pool
import logging
import string
import os, os.path

__program__ = 'run_ivy'
__author__ = 'Soh Ishiguro <yukke@g-language.org>'
__license__ = ''
__status__ = 'development'
__version__ = __version__

def run():
    logger = logging.getLogger(__name__)    
    parse = CommandLineParser()
    params = parse.ivy_parse_options()
    print params
    vcf = VCFWriteHeader(params)
    #vcf.make_vcf_header()
    logger.debug("Beginning Ivy run (v." + __version__ + ")" )
    
    if params.r_bams:
        logger.debug("Loading RNA-seq bam file '{0}'".format(params.r_bams))
        rna_pileup_alignment = RNASeqAlignmentStream(params)
        for rna in rna_pileup_alignment.filter_stream():
            print to_tab(rna)
            
    if params.d_bams:
        logger.debug("Loading DNA-seq bam file '{0}'".format(params.d_bams))
        dna_pileup_alignment = DNASeqAlignmentStream(params)
        for dna in dna_pileup_alignment.filter_stream():
            print '{chrom:}\t{pos:}\t{ref:}\t{alt:}'.format(
                chrom=dna.get('chrom'),
                pos=dna.get('pos'),
                ref=dna.get('ref'),
                alt=dna.get('alt'))
    #with open(params.outname, 'w') as f:
    #f.write()

def __thread_ivy():
    parse = CommandLineParser()
    params = parse.ivy_parse_options()
    vcf = VCFWriteHeader(params)
    if params.r_bams:
        pileup_iter = RNASeqAlignmentStream(params)
        for pileup in pileup_iter.filter_stream():
            print to_tab(pileup)
    elif params.d_bams:
        pileup_iter = DNASeqAlignmentStream(params)
        for pileup in pileup_iter.filter_stream():
            print to_tab(pileup)
    
def __worker(cpus=1, seqs=None):
    if cpus <= 1 and len(seqs) <= 1:
        raise RuntimeError()
    p = Pool(cpus)
    seq = p.map(__thread_ivy, seqs)
    return seq

def thread_run():
    parse = CommandLineParser()
    params = parse.ivy_parse_options()
    path = './block_fasta'
    fasta_file = params.fasta
    cpus = params.n_threads
    
    # create tmp directory
    if not os.path.isdir(path):
        os.mkdir(path)
    fasta_files = [_ for _ in os.listdir(path) if _.endswith('.fa')]
    
    # generate worker
    if len(fasta_files) != cpus:
        print "Remove old splited fasta, and generate new..."
        shutil.rmtree(path)
        fa = Fasta(fa=fasta_file)
        fa.split_by_blocks(n=cpus)
        __worker(cpus=cpus, seqs=fasta_files)
        
    elif len(fasta_files):
        print "Used existing splited fasta..."
        __worker(cpus=cpus, seqs=fasta_files)
        
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
        
    
