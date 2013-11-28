import string
import re
import os.path
import multiprocessing
import pprint

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

    def chr_size(self):
        return len(self.fasta_header())
    
    def split_by_length(self, length):
        for i in range(1, num):
            f = open(str(i)+'_'+str(self.filename), 'w')
            f.write(str(i))
            f.close()
        
    def split_by_cpus(self, cpus):
        debug = True
        
        #MAX_CPUs = multiprocessing.cpu_count()
        MAX_CPUs = 24
        if cpus > MAX_CPUs:
            raise RuntimeError("Over the number of cpus are given")
            
        human_chr = ['chrM', 'chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6',
                     'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13',
                     'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20',
                     'chr21', 'chr22', 'chrX', 'chrY']
        
        num_threads = 9
        chr_size = len(human_chr)
        div, mod = divmod(chr_size, num_threads)
        start = 0
        end = div
            
        result = []
        overflow = []
        
        if debug:
            print "Num threads: ", num_threads
            print "div: ", div
            print "mod: ", mod
            
        counter = num_threads
        for i in range(0, chr_size):
            if len(human_chr[start:end]):
                if mod == 0:
                    result.append(human_chr[start:end])
                    end += div
                    start += div
                elif mod >= 1:
                    if counter > 0:
                        # divisible
                        result.append(human_chr[start:end])
                        start += div
                        end += div
                    elif counter == 0:
                        # has mod
                        overflow.append(human_chr[start:end])
                counter -= 1
        
        def __merge_list(norm, over):
            for index, elem in enumerate(norm):
                if index < len(over[0]):
                    norm[index].append(over[0][index])
            return norm
        return __merge_list(result, overflow)

        
if __name__ == '__main__':
    fasta = Fasta(fa="/Users/yukke/dev/data/genome.fa")
    #print fasta.fasta_header()
    #print fasta.split_by_chr)
    
    print fasta.split_by_cpus(5)
