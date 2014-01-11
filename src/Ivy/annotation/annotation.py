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
    def __init__(self, ingtf, contig=None, start=None, end=None):
        self.contig = contig
        self.start = start
        self.end = end
        self._prepare()
        
    def _prepare(self):
        if not os.path.isfile(ingtf + ".gz.tbi"):
            print "Indexing GTF file..."
            compressed_gtf = pysam.tabix_index(ingtf, preset="gff")
        else:
            compressed_gtf = ingtf + ".gz"
        self.tabixfile = pysam.Tabixfile(compressed_gtf)

    def fetch_gtf(self):
        for gtf in pysam.Tabixfile.fetch(self.tabixfile,
                                         parser=pysam.asGTF(),
                                         start=self.start, end=self.end):
            yield gtf
            
if __name__ == '__main__':
    pass
