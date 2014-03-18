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

    def gene_in_region(self, contig=None, start=None, end=None):
        for gtf in pysam.Tabixfile.fetch(self.tabixfile, contig, start, end,
                                         parser=pysam.asGTF()):
            yield gtf.asDict()['gene_name']

    def subset_of_feature_in_region(self, contig=None, start=None, end=None, types=None):
        '''
        Example:
         # return dict of rna and tss id
         for rna_tss in tabix.subset_of_feature_in_region(contig="chr2L", end=9839,
                                                    types=["transcript_id", 'tss_id']):
             print rna_tss
        '''
        for gtf in pysam.Tabixfile.fetch(self.tabixfile, contig, start, end,
                                         parser=pysam.asGTF()):
            if isinstance(types, str):
                try:
                    yield gtf.asDict()[types]
                except KeyError:
                    print 'key \'{0}\' is not found in {1}'.format(t, self.ingtf)
                    
            elif isinstance(types, list):
                tmp = dict()
                for t in types:
                    try:
                        tmp.update({t: gtf.asDict()[t]})
                    except KeyError:
                        print 'key \'{0}\' is not found in {1}'.format(t, self.ingtf)
                yield tmp

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
    for _ in tabix.subset_of_feature_in_region(contig="chr2L", end=9839, types=["transcript_id", 'tss_id']):
        print _
        #print _.asDict()
        #print  dir(_)
        #print _.asDict()
        
