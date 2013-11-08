#!/usr/bin/env python

from Ivy.alignment.stream import *
from Ivy.annotation.writer import VCFWriteHeader
import sys
import pprint
import csv

if __name__ == '__main__':
    SMALL = False
    
    if SMALL:
        bam_file = '/Users/yukke/dev/data/testREDItools/rna.bam'
        fa_file = '/Users/yukke/dev/data/testREDItools/reference.fa'
        chr_name = 'chr21'
        start = 47721030
        end = 47721057
        
    else:
        bam_file = '/Users/yukke/dev/data/GSM958732_hg19_wgEncodeCaltechRnaSeqHepg2R2x75Il200SplicesRep1V2.sorted.bam'
        fa_file = '/Users/yukke/dev/data/testREDItools/reference.fa'
        #fa_file = '/Users/yukke/dev/data/genome.fa'
        chr_name = 'chr21'
        start = 47721030
        end =   48129895
        
    rna_alignment = AlignmentStream(bam_file, fa_file, chrom=chr_name, start=start, end=end)
    
    pp = pprint.PrettyPrinter(indent=6)
    
    for rna in rna_alignment.pileup_stream():
        #print "{0}\t{1}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\n".format(rna['CHROM'], rna["POS"], rna["ID"], rna["REF"], rna["ALT"], rna["QUAL"], rna["mismatch_freq"], rna["coverage"], rna["FORMAT"])

        #print "%s\t%s\t%s\t%s\t%s\t%s\t%s" % \
            #(rna['CHROM'], rna['POS'], rna['REF'], rna['raw_coverage'], rna['nodel_coverage'], rna['prop_coverage'], rna['prop_nodel_coverage'])
        print [str(_)+ ":"+ str(rna[_])for _ in rna]
