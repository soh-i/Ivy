from multiprocessing import Pool, Process

chroms = ['chr'+str(_) for _ in range(1, 23)]
chroms += ['chrX', 'chrY']

threads = 7
max_cpus = 24

try:
    div, mod = divmod(len(chroms), threads)
except ZeroDivisionError:
    pass

print "threads:", threads
print div, mod

c = 1
start = 0
end = div
result = []
overflow = []
for num in range(1, len(chroms)+1):
    if threads > max_cpus:
        raise RuntimeError("Threads num is bigger than system CUPs")
    mm = []
    if len(chroms[start:end]) > 0:
        if div != 0 and mod == 0:
            #print num, chroms[start:end]
            start += div
            end += div
        elif div != 0 and mod != 0:
            if c < (len(chroms) - mod) and num < threads:
                result.append(chroms[start:end])
                #print num, chroms[start:end]
                start += div
                end += div
            elif c == threads:
                # over flow
                result.append(chroms[start:end])
                #print  num, chroms[start:end]
                start += div
                end += div
                overflow.append(chroms[start:end])
            c += 1

for i in range(0, mod):
    result[i].extend([overflow[0][i]])

print result


"""
def split_fasta(fastafile, num):
    if not os.path.isfile(fastafile):
        raise ValueError("{0} is not found".format(fastafile))
        
            

class test(object):
    def __init__(self, x):
        self.x = x

    def time(self):
        return self.x ** self.x

def run(x):
    t = test(x)
    t.time()
    
if __name__ == '__main__':
    p = Pool(4)
    #op = p.map(run, range(10000))
    p = Process(target=run, args=(100000000,))
    p.start()
    p.join()
    
    

"""
