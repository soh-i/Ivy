from Ivy.benchmark import *

if __name__ == '__main__':
    darned_reader = DarnedReader(sp='human')
    v = VCFReader('/home/soh.i/RNA_editing/sandbox/rec3h_205_reed0_mis2gap0dist0.vcf')
    pred = v.generate_vcf_set()
    ans = darned_reader.generate_darned_set()
    bench = Benchmark(answer=ans, predict=pred)
    print bench.precision()
    print bench.recall()
    print bench.f_measure()

    

    
    
    
