#!/usr/bin/env python

from Ivy.alignment import *

if __name__ == '__main__':
    bam_file = '../test_data/dna.bam'
    #bam_file = '/home/soh.i/db/melanogaster/Nascent-Seq/ZT18_R1/accepted_hits.bam'
    fa_file = '../test_data/reference.fa'
    #fa_file = '/home/soh.i/archives/nucRNA/data/Drosophila_melanogaster/UCSC/dm3/Sequence/WholeGenomeFasta/genome.fa'
    chr_name = 'chr21'
    start = 0
    end = 0
    alignment = Alignment(bam_file, fa_file, chrom=chr_name, start=start, end=end)
    
    for i in alignment.pileup_stream():
        print i['chrom'], i['pos'], i['ref'],
        print i['coverage'],
        print i['mismatches'],
        print i['Af'],
        print i['Ar'],
        print i['Tf'],
        print i['Tr'],
        print i['Gf'],
        print i['Gr'],
        print i['Cf'],
        print i['Gr'],
        print i['N']
        
