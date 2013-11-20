from Ivy.version import __version__
from Ivy.utils import die
from Ivy.alignment.stream import AlignmentStream
from Ivy.parse_opt import CommandLineParser
from Ivy.annotation.writer import VCFWriteHeader
import pprint

__program__ = 'run_ivy'
__author__ = 'Soh Ishiguro <yukke@g-language.org>'
__license__ = ''
__status__ = 'development'


class Ivy(object):
    def __init__(self):
        pass
        
    def run(self):
        pp = pprint.PrettyPrinter(indent=6)
        parse = CommandLineParser()
        params = parse.ivy_parse_options()
        vcf = VCFWriteHeader(params)
        stream = AlignmentStream(params)
        vcf.make_vcf_header()
        
        for rna in stream.pileup_stream():
            print '{chrom}\t{pos:d}\t{ref:s}\t{alt:s}\t{coverage:d}\t{mismatch_ratio:0.4f}\t{dp4:s}'.format(
                chrom=rna['chrom'],
                pos=rna['pos'],
                ref=rna['ref'],
                alt=rna['alt'],
                coverage=rna['coverage'],
                mismatch_ratio=rna['mismatch_ratio'],
                dp4=":".join([str(_) for _ in rna['dp4']]))
