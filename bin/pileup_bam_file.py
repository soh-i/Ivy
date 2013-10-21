import vcf
import pysam

def mismatch_stat(samfile, fafile, chrom=None, start=None, end=None, one_based=True):
    coords = resolve_coords(one_based, start, end)
    
    for col in samfile.pileup(reference=chrom, start=coords['start'], end=coords['end']):
        chrom = samfile.getrname(col.tid)
        pos = col.pos + 1 if one_based else col.pos

        prop_read = []
        dele_read = []
        inser_read = []
        
        for r in col.pileups:
            if r.alignment.is_proper_pair and not r.alignment.is_duplicate and not r.alignment.is_unmapped and not r.is_del:
                prop_read.append(r)
                
        ref = fafile.fetch(chrom, col.pos, col.pos+1).upper()
        mismatches = [read for read in prop_read
                      if read.alignment.seq[read.qpos] != ref]
        matches = [read for read in prop_read
                   if read.alignment.seq[read.qpos] == ref]
        
        A = [read for read in prop_read
             if read.alignment.seq[read.qpos] == 'A']
        C = [read for read in prop_read
             if read.alignment.seq[read.qpos] == 'C']
        T = [read for read in prop_read
             if read.alignment.seq[read.qpos] == 'T']
        G = [read for read in prop_read
             if read.alignment.seq[read.qpos] == 'G']
        N = [read for read in prop_read
             if read.alignment.seq[read.qpos] == 'N']

        yield {'chrom': chrom,
               'pos': pos,
               'ref': ref,
               'proper_reads': len(prop_read),
               'matches': len(matches),
               'mismatches': len(mismatches),
               'A': len(A),
               'C': len(C),
               'T': len(T),
               'G': len(G),
               'N': len(N)
               }

def resolve_coords(one_based, start, end):
    coords = {'start':None, 'end':None}
    if one_based:
        if start is not None:
            coords.update({'start':start-1})
        if end is not None:
            coords.update({'end':end-1})
    return coords

if __name__ == '__main__':
    bam_file = '../test_data/dna.bam'
    bm = pysam.Samfile(bam_file, 'rb')
    fa_file = '../test_data/reference.fa'
    ft = pysam.Fastafile(fa_file)
    chr_name = 'chr21'
    start = 47743799
    end =   47743800

    for i in mismatch_stat(bm, ft, chrom=chr_name, start=start, end=end):
        print i['chrom'], i['pos'], i['ref'],
        print i['matches'],
        print i['A'],
        print i['T'],
        print i['G'],
        print i['C']
        
