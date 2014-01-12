import pysam
import os.path

class Annotation(object):
    def to_array(self):
        pass

    def to_file(self):
        pass
        

class AnnotateVCF(Annotation):
    def __init__(self, vcf=None, gtf=None):
        self.__vcf = vcf
        #self.__gtf = gtf

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
                    

    #def gtf_to_array(self):
    #    with open(self.__gtf, 'w') as f:
    #        for line in f:
    #            data = split.line("\t")
                
    def compare(self):
        pass

    def write_gtf(self):
        pass


class GTF(Annotation):
    def __init__(self, ingtf):
        self.ingtf = ingtf
        self._prepare()
        
    def _prepare(self):
        if not os.path.isfile(self.ingtf + ".gz.tbi"):
            print "Indexing GTF file..."
            compressed_gtf = pysam.tabix_index(self.ingtf, preset="gff")
        else:
            compressed_gtf = self.ingtf + ".gz"
        self.tabixfile = pysam.Tabixfile(compressed_gtf)

    def fetch_gtf(self, contig=None, start=None, end=None):
        for gtf in pysam.Tabixfile.fetch(self.tabixfile, contig, start, end,
                                         parser=pysam.asGTF()):
                                         
            yield gtf

    def strand_info(self, contig=None, start=None, end=None):
        for gtf in pysam.Tabixfile.fetch(self.tabixfile, contig, start, end,
                                         parser=pysam.asGTF()):
            return gtf.strand

if __name__ == '__main__':
    gtf_cls = GTF("/Users/yukke/Desktop/genes.gtf")
    print gtf_cls.strand_info("chr21", 48084236, 49094239)
    
        

