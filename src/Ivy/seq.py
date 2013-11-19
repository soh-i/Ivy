import string
import re
import os.path

class Fasta(object):
    def __init__(self, fa=None):
        if os.path.isfile(fa):
            self.filename = fa
        else:
            raise RuntimeError("{fa:s} is not found".format(fa=fastafile))
    
    def fasta_header(self):
        flag = False
        count = 0
        header = []
        with open(self.filename, 'r') as f:
            for line in f:
                if line.startswith('>'):
                    header.append(line.rstrip())
                    count += 1
                    flag = True
        return header

    def split_fasta(self, num):
        for i in range(1, num):
            f = open(str(i)+'_'+str(self.filename), 'w')
            f.write(str(i))
            f.close()
        
    def get_fasta_length(self):
        raise NotImplementedError()
        

if __name__ == '__main__':
    fa = Fasta(fa="seq.fa")
    print fa.filename
    print fa.fasta_header()
    print fa.split_fasta(4)
