from Ivy.benchmark import *

if __name__ == '__main__':
    darned_reader = DarnedReader(sp='human_hg19')
    print darned_reader
    v = VCFReader('../data/test_data.vcf')
    pred = v.generate_vcf_set()
    ans = darned_reader.generate_darned_set()
    bench = Benchmark(answer=ans, predict=pred)

    
    #print "Precision:\t%f" % (bench.precision())
    #print "Recall:\t%f" % (bench.recall())
    #print "F-measure:\t%f" % (bench.f_measure())

