#!/usr/bin/env python

from Ivy.alignment import *
from Ivy.writer import VCFWriteHeader

if __name__ == '__main__':
    bam_file = '../test_data/dna.bam'
    #bam_file = '/home/soh.i/db/melanogaster/Nascent-Seq/ZT18_R1/accepted_hits.bam'
    fa_file = '../test_data/reference.fa'
    #fa_file = '/home/soh.i/archives/nucRNA/data/Drosophila_melanogaster/UCSC/dm3/Sequence/WholeGenomeFasta/genome.fa'
    chr_name = 'chr21'
    start = 47721045
    end = 47721088
    
    rna_alignment = AlignmentStream(bam_file, fa_file, chrom=chr_name, start=start, end=end)
    #dna_alignment = AlignmentStream(bam_file, fa_file, chrom=chr_name, start=start, end=end)

    v = VCFWriteHeader()
    v.make_vcf_header()
   
    for rna in rna_alignment.pileup_stream():
        print rna

"""        
        for dna in dna_alignment.pileup_stream():
            if rna['chrom'] == dna['chrom'] and dna['pos'] == rna['pos']:
                print dna['pos']

                break;
                
            
"""
