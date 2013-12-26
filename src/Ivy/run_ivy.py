from Ivy.version import __version__
from Ivy.utils import die
from Ivy.alignment.stream import AlignmentStream, RNASeqAlignmentStream, DNASeqAlignmentStream
from Ivy.commandline.parse_ivy_opts import CommandLineParser
from Ivy.annotation.writer import VCFWriteHeader
from  pprint import pprint as p

__program__ = 'run_ivy'
__author__ = 'Soh Ishiguro <yukke@g-language.org>'
__license__ = ''
__status__ = 'development'
__version__ = __version__

def run():
    parse = CommandLineParser()
    params = parse.ivy_parse_options()
    vcf = VCFWriteHeader(params)
    #stream = AlignmentStream(params)
    rna_stream = RNASeqAlignmentStream(params)
    dna_stream = DNASeqAlignmentStream(params)
    vcf.make_vcf_header()

    #with open(params.outname, 'w') as f:
    #f.write()
            
    for rna in rna_stream.pileup_stream():
        pprint(rna)
        #print '{chrom}\t{pos:d}\t{ref:s}\t{alt:s}\tDP={coverage:d};DP4={dp4:s};MIS_RATIO={mismatch_ratio:0.4f}'.format(
        #    chrom=rna['chrom'],
        #    pos=rna['pos'],
        #    ref=rna['ref'],
        #    alt=rna['alt'],
        #    mismatch_ratio=rna['mismatch_ratio'],
        #    coverage=rna['coverage'],
        #    dp4=",".join([str(_) for _ in rna['dp4']]))
        
    for dna in dna_stream.pileup_stream():
        pprint(dna)


def pprint(data, *args, **kwargs):
    '''
    Simple pretty print yielded pileuped data
    '''
    
    print '{'
    for key in data:
        print '  {key}({types}): {val}'.format(
            key=key, val=data[key], types=type(data[key]).__name__)
    print '}'
