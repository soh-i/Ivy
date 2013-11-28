import string
import re
import os.path
import multiprocessing

class Fasta(object):
    def __init__(self, fa=''):
        if os.path.isfile(fa):
            self.filename = fa
        else:
            raise RuntimeError("{:s} is not found".format(fa))
    
    def fasta_header(self):
        header = []
        with open(self.filename, 'r') as f:
            #head = [line.replace('>', '', 1).rstrip() for line in f if line.startswith('>')]
            for line in f:
                if line.startswith('>'):
                    head = line.replace('>', '', 1).rstrip()
                    header.append(head)
            return header

    def split_by_chr(self):
        flag = False
        count = 0
        with open(self.filename, 'r') as f:
            for line in f:
                if line.startswith('>'):
                    count += 1
                    head = line.replace('>', '', 1).rstrip()
                    out = open(head + '.fa', 'w')
                    flag = True
                if flag:
                    out.write(line)
                elif not flag:
                    continue
        return None
                
    def split_by_length(self, length):
        for i in range(1, num):
            f = open(str(i)+'_'+str(self.filename), 'w')
            f.write(str(i))
            f.close()
        
    def get_fasta_length(self):
        pass

def automatically(cpus):
    #MAX_CPUs = multiprocessing.cpu_count()
    MAX_CPUs = 24
    if cpus > MAX_CPUs:
        raise RuntimeError("Over the number of cpus are given")
        

if __name__ == '__main__':
    fasta = Fasta(fa="/Users/yukke/dev/data/genome.fa")
    #print fasta.fasta_header()
    #print fasta.split_by_chr()

    
    print automatically(24)

    
