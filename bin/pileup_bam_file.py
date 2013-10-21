import vcf
import pysam

bam_file = '../test_data/dna.bam'
bam = pysam.Samfile(bam_file, 'rb')
fa_file = '../test_data/reference.fa'
fasta = pysam.Fastfile(fa_file)

chr_name = 21
start = 47728970
end = 47745353

for column in bam.pileup(chr_name, start, end):
    print column
    

