#!/usr/bin/env python

from Ivy.alignment import *
from Ivy.writer import VCFWriteHeader
import sys
import pprint

if __name__ == '__main__':
    SMALL = not True
    bam_file = fa_file = chr_name = start = end = 0
    
    if SMALL:
        bam_file = 'test_data/dna.bam'
        fa_file = 'test_data/reference.fa'
        chr_name = 'chr21'
        start = 1
        end = 100
        
    else:
        bam_file = '/home/soh.i/db/melanogaster/Nascent-Seq/ZT18_R1/accepted_hits.bam'
        fa_file = '/home/soh.i/archives/nucRNA/data/Drosophila_melanogaster/UCSC/dm3/Sequence/WholeGenomeFasta/genome.fa'
        chr_name = '2L'
        start = 1
        end = 1000

    ft = pysam.Fastafile(fa_file)
    ref = ft.fetch(reference=chr_name, start=start, end=end)
    
    if ref:
        print True
    else:
        print False
    
    rna_alignment = AlignmentStream(bam_file, fa_file, chrom=chr_name, start=start, end=end)
    #dna_alignment = AlignmentStream(bam_file, fa_file, chrom=chr_name, start=start, end=end)

    v = VCFWriteHeader()
    v.make_vcf_header()
    pp = pprint.PrettyPrinter(indent=6)
    
    for rna in rna_alignment.pileup_stream():
        #print "{0}\t{1}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\n".format(rna['CHROM'], rna["POS"], rna["ID"], rna["REF"], rna["ALT"], rna["QUAL"], rna["mismatches"], rna["coverage"], rna["FORMAT"]),
        pp.pprint(rna)
        

