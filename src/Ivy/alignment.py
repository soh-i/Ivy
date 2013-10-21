import pysam

class Alignment(object):
    def __init__(self, samfile, fafile, chrom=None, start=None, end=None, one_based=None):
        bm = pysam.Samfile(samfile, 'rb')
        ft = pysam.Fastafile(fafile)
        
        self.samfile = bm
        self.fafile = ft
        self.chrom = chrom
        self.start = start
        self.end = end
        self.one_based = one_based
        
    def to_pileup(self):
        coords = self.__resolve_coords()
    
        for col in self.samfile.pileup(reference=self.chrom,
                                       start=self.coords['start'],
                                       end=self.coords['end']):
            chrom = self.samfile.getrname(col.tid)
            pos = col.pos + 1 if self.one_based else col.pos

            prop_read = []
            for r in col.pileups:
                if r.alignment.is_proper_pair \
                   and not r.alignment.is_duplicate \
                   and not r.alignment.is_unmapped \
                   and not r.is_del:
                    prop_read.append(r)
                
            ref = self.fafile.fetch(self.chrom, col.pos, col.pos+1)
            mismatches = [read for read in prop_read
                          if read.alignment.seq[read.qpos] != ref]
            matches = [read for read in prop_read
                       if read.alignment.seq[read.qpos] == ref]
        
            A = [read for read in prop_read
                 if read.alignment.seq[read.qpos] == 'A']
            a = [read for read in prop_read
                 if read.alignment.seq[read.qpos] == 'n']
            
            C = [read for read in prop_read
                 if read.alignment.seq[read.qpos] == 'C']
            c = [read for read in prop_read
                 if read.alignment.seq[read.qpos] == 'c']
            
            T = [read for read in prop_read
                 if read.alignment.seq[read.qpos] == 'T']
            t = [read for read in prop_read
                 if read.alignment.seq[read.qpos] == 't']
            
            G = [read for read in prop_read
                 if read.alignment.seq[read.qpos] == 'G']
            g = [read for read in prop_read
                 if read.alignment.seq[read.qpos] == 'g']
            
            N = [read for read in prop_read
                 if read.alignment.seq[read.qpos] == 'N' \
                 or read.alignment.seq[read.qpos] == 'n']

            yield {'chrom':chrom,
                   'pos':pos,
                   'ref':ref,
                   'coverage':len(prop_read),
                   'matches':len(matches),
                   'mismatches':len(mismatches),
                   'Af':len(A),
                   'Ar':len(a),
                   'Cf':len(C),
                   'Cr':len(c),
                   'Tf':len(T),
                   'Tr':len(t),
                   'Gf':len(G),
                   'Gr':len(g),
                   'N': len(N)
               }

    def __resolve_coords(self):
        self.coords = {'start':None, 'end':None}
        if self.one_based:
            if self.start is not None:
                self.coords.update({'start':self.start-1})
            if end is not None:
                self.coords.update({'end':self.end-1})
        return self.coords

