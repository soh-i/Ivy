import pysam

class Annotation(object):
    def __init__(self):
        self.gtf = gtf

    def to_array(self):
        pass
        
class AnnotateVCF(object):
    def __init__(self, vcf=None, gtf=None):
        self.__vcf = vcf
        #self.__gtf = gtf

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

        
def main():
    a = AnnotateVCF(vcf="../../../data/test_data.vcf")

    print a.vcf_to_array()[:10][:10]

if __name__ == '__main__':
    main()




    
