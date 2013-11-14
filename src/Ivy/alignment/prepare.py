import os.path
import pysam

__program__ = 'prepare'
__author__ = 'Soh Ishiguro <yukke@g-language.org>'
__license__ = ''
__status__ = 'development'


class AlignmentPreparation(object):
    def __init__(self):
        # TODO: AlignmentStream.__init__ move into here.
        pass
        
    def alignment_prepare(self):
        raise NotImplementedError
    
    def __sort(self):
        if not os.path.isfile(bamfile):
            try:
                pysam.sort(self.samfile, self.samfile + 'sorted')
                sort_log = pysam.sort.getMessage()
                return True
            except:
                raise RuntimeError()
        else:
            print "already sorted"
            return False

    def __index(self):
        if not os.path.isfile(samfile + '.index.bam'):
            try:
                pysam.index(self.samfile)
                return True
            except:
                raise RuntimeError()
        else:
            print "already indexed"
            return False

    def __faidx(self):
        if not os.path.isfile(fafile + '.fai'):
            try:
                pysam.faidx(self.fafile)
                return True
            except:
                raise RuntimeError()
        else:
            print "already exist"
            return False 

    def __merge_bams(self, bams=[]):
        for _ in bams:
            if not os.path.isfile(_):
                raise RuntimeError()
        try:
            pysam.merge([_ for _ in bams])
            return True
        except:
            raise RuntimeError()
