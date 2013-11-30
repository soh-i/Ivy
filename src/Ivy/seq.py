import string
import re
import string
import os.path
import os
import shutil
import multiprocessing
import pysam
import pprint
import time

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

    def split_by_blocks(self, n=0):
        '''
        Args:
         lists(list): block of splited chromosome
        '''
        
        self.save_path = './block_fasta/'
        if not os.path.isdir(self.save_path):
            os.mkdir(self.save_path)
            
        block_list = self.generate_chrom_blocks(n)
        count = 0
        block_size = len(block_list)
        if block_size < 0:
            raise ValueError("Block size must be greater than 1")
        
        for block in block_list:
            for chrom in block:
                name = "-".join(block)
                out = open(self.save_path + str(block_size) + '_' + name + '.fa', 'w')
                fa = open(self.filename, 'r')
                
                for line in fa:
                    if line.startswith('>'):
                        head = line.replace('>', '', 1).rstrip()
                    if head in block:
                        out.write(line)
                break
                out.close()
                fa.close()
            block_size -= 1

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
            
        #human_chr = ['chrM', 'chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6',
        #             'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13',
        #             'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20',
        #             'chr21', 'chr22', 'chrX', 'chrY']
        
        human_chr = self.fasta_header()
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


def decode_chr_name_from_file(chroms):
    '''
    Args:
     blocked fasta file, 1_chr18-chr19-chr20.fa
    Returns:
     chromosome name(list): [chr18, chr19, chr20]
    '''
    decode = []
    files = [f.split("-") for f in chroms]
    for i in files:
        for chrm in i:
            # remove prefix
            p = re.compile(r'(^\d+_)(.+)')
            match = p.match(chrm)
            if match:
                decode.append(match.group(0).replace(match.group(1), '', 1))
            else:
                # remove suffix
                p = re.compile(r'(.+)\.{1}fa$')
                match = p.match(chrm)
                if match:
                    decode.append(match.group(1))
                else:
                    # normal
                    decode.append(chrm)
    return decode
                
def fetch_seq(fa):
    '''
    Args:
     fasta(list): Path to fasta file
    '''
    
    print multiprocessing.current_process()
    
    path = './block_fasta/'
    fafile = pysam.Fastafile(path+fa)
    chroms = decode_chr_name_from_file([fa])
    seq = []
    for c in chroms:
        #print "fetching {0}...".format(c)
        l = len(fafile.fetch(reference=c, start=1, end=1000000000))
        print "Chr: %s, Length: %d" % (c, l)
        
    return seq

def run(cpus, fas):
    '''
    Args:
     cpus(int): Number of cpus
     fasta(list): Path to fasta files
    '''
    
    p = multiprocessing.Pool(cpus)
    seq = p.map(fetch_seq, fas)
    return seq

def get_fa_list(path):
    '''
    Args:
     path(str): blocked fasta contained directory

    Returns:
     files(list): only .fa file which located in given path 
    '''

    fa = []
    for _ in os.listdir(path):
        if _.endswith(".fa"):
            fa.append(_)
    return fa

class Timer(object):
    def __init__(self, verbose=True):
        self.verbose = verbose
        
    def __enter__(self):
        self.start = time.time()

    def __exit__(self, *args):
        self.end = time.time()
        self.secs = self.end - self.start
        if self.verbose:
            print '=> Elapsed time: {:f} (sec)'.format(self.secs)
    
if __name__ == '__main__':
    path = './block_fasta/'
    fasta_file = "../../../data/genome.fa"
    cpus = 4
    
    # create directory
    if not os.path.isdir(path):
        os.mkdir(path)
    fas = get_fa_list(path)
    if len(fas) != cpus:
        print "Remove old blocks and generates new blocks"
        shutil.rmtree(path)
        fa = Fasta(fa=fasta_file)
        fa.split_by_blocks(n=cpus)
        
        with Timer() as t:
            run(cpus, get_fa_list(path))
        
    elif len(get_fa_list(path)) == cpus:
        print "Used existing block"
        with Timer() as t:
            run(cpus, get_fa_list(path))
            
    ##start = time.clock()
    #for f in fas:
    #    fafile = pysam.Fastafile(path+f)
    #    seq = fafile.fetch(reference=decode_chr_name(f),start=1, end=100000)
    #    print fafile.filename
    #    #print seq
    #end = time.clock()
    #print end - start
    #    
