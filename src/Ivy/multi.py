from multiprocessing import Pool, Process

chroms = ['chr'+str(_) for _ in range(1, 23)]
chroms += ['chrX', 'chrY']

threads = 5
max_cpus = 24
div, mod = divmod(len(chroms), threads)
start = 0
end = div
result = []
print "threads:", threads
print div, mod

c = 1
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
                mm.append(chroms[start:])
                
                #print num, chroms[start:]
            c += 1


    for _ in mm:
        if _:
            for i in _:
                #result.append(i)
                print i[3]
#print resul
            
def split_fasta(fastafile, num):
    if not os.path.isfile(fastafile):
        raise ValueError("{0} is not found".format(fastafile))
        
            
"""
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
