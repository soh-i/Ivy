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

    def split_by_blocks(self, lists):
        flag = False
        count = 0
        print lists
        
        with open(self.filename, 'r') as f:
            for line in f:
                if line.startswith('>'):
                    count += 1
                    head = line.replace('>', '', 1).rstrip()
                    out = open(str(count) + ".fa", 'w')
                    flag = True
                if flag:
                    #out.write(line)
                    print line,
                elif not flag:
                    continue
        return None

    def chr_size(self):
        return len(self.fasta_header())
    
    def split_by_length(self, length):
        pass

    def generate_chrom_blocks(self, cpus):
        '''
        Args:
         cpus(int): number of cpus

        Returns:
         chrom(list): splited chromsome by cpus
        '''
        
        #MAX_CPUs = multiprocessing.cpu_count()
        MAX_CPUs = 24
        if cpus > MAX_CPUs:
            raise RuntimeError("Over the number of cpus are given")
            
        human_chr = ['chrM', 'chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6',
                     'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13',
                     'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20',
                     'chr21', 'chr22', 'chrX', 'chrY']
        
        num_threads = cpus
        chr_size = len(human_chr)
        try:
            div, mod = divmod(chr_size, num_threads)
        except ZeroDivisionError as e:
            raise SystemExit(e)
            
        start = 0
        end = div
        result, overflow = [], []
        counter = num_threads
        
        for i in range(0, chr_size):
            if len(human_chr[start:end]):
                if mod == 0:
                    result.append(human_chr[start:end])
                    end += div
                    start += div
                    counter -= 1
                elif mod >= 1:
                    if counter > 0:
                        result.append(human_chr[start:end])
                        start += div
                        end += div
                        counter -= 1
                    else:
                        overflow.append(human_chr[start:end])
                        start += div
                        end += div
                        counter -= 1
                        
        def __merge_list(norm, over):
            flatten = reduce(lambda x, y: x + y, over)
            for i, _ in enumerate(norm):
                for j, _ in enumerate(flatten):
                    norm[i].append(flatten[j])
                    del flatten[j]
                    break
            return norm

        if len(overflow) > 0:
            return __merge_list(result, overflow)
        elif len(overflow) == 0:
            return result

        
if __name__ == '__main__':
    #fasta = Fasta(fa="/Users/yukke/dev/data/genome.fa")
    fasta = Fasta(fa="seq.fa")
    blocks = fasta.generate_chrom_blocks(5)
    fasta.split_by_blocks(blocks)
    

    
