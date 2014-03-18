import pysam
import os.path
import logging

class AnnotateVCF(object):
    def __init__(self, vcf=None, gtf=None):
        self.__vcf = vcf
        self.nuc = {'A': 'T',
                    'T': 'A',
                    'C': 'G',
                    'G': 'C',
                    }
        
    def vcf_to_array(self):
        recs = []
        with open(self.__vcf, 'r') as f:
            for line in f:
                if not line.startswith("#"):
                    data = line.split("\t")
                    chrom = data[0]
                    pos = data[1]
                    recs.append([chrom, pos])
        return recs
                    

class GTF(object):
    def __init__(self, ingtf):
        self.ingtf = ingtf
        self._prepare()
        
    def _prepare(self):
        if not os.path.isfile(self.ingtf + ".gz.tbi"):
            print "Generate indexed GTF (tabix) file: '{0}'...".format(self.ingtf)
            compressed_gtf = pysam.tabix_index(self.ingtf, preset="gff")
        else:
            compressed_gtf = self.ingtf + ".gz"
        self.tabixfile = pysam.Tabixfile(compressed_gtf)

    def fetch_gtf(self, contig=None, start=None, end=None):
        '''
        Returns:
         pysam.tabix object
        
        Example:
         # return dict of each GTF line
         for _ in tabix.fetch_gtf(contig="chr2L", end=9839):
             print _.asDict()
        '''
        
        for gtf in pysam.Tabixfile.fetch(self.tabixfile, contig, start, end,
                                         parser=pysam.asGTF()):
            yield gtf

    def strand_info(self, contig=None, start=None, end=None):
        """
        Args:
         contig(str)='', start(int)='', end=''
        
        Returns:
         strand information [+-], or [.] is 404
        """
        
        debug = False
        found = "."
        for gtf in pysam.Tabixfile.fetch(self.tabixfile, contig, start, end,
                                         parser=pysam.asGTF()):
            if gtf.strand:
                found = gtf.strand
                if debug:
                    print gtf.gene
                    print "start: %s, end: %s" % (gtf.start, gtf.end)
                break
            else:
                continue
        return found

if __name__ == '__main__':
    tabix = GTF("genes.gtf")
    for _ in tabix.fetch_gtf(contig="chr2L", end=9839):
        print _.asDict()
        
