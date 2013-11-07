from Ivy.benchmark.benchmark import (
    DarnedDataGenerator,
    DarnedReader,
    VCFReader,
    Benchmark,
    )
from Ivy.benchmark.plot import BenchmarkPlot
import Ivy.utils
import argparse
import os.path


def run():
    args = parser.parse_args()
    gen = DarnedDataGenerator(species=args.sp)
    
    darned_raw_file = gen.saved_abs_path + gen.filename
    if not os.path.isfile(darned_raw_file):
        print "fetching from darned..."
        gen.fetch_darned()
        
    darned_parsed_csv = gen.out_name
    if not os.path.isfile(darned_parsed_csv):
        print "parsing darned..."
        gen.darned_to_csv()
        
    if args.vcf_file and args.sp:
        ans = DarnedReader(sp=args.sp, source=args.source)
        vcf = VCFReader(args.vcf_file)
        bench = Benchmark(answer=ans.db, predict=vcf.db)
        
        print "Species:%s,DB:%s,VCF:%s,Precision:%f,Recall:%f,F-measure:%f,AGs:%d,Others:%d,AnsCount:%d" % (
            ans.sp()[0], ans.db_name(), vcf.vcf_name(),
            bench.precision(), bench.recall(), bench.f_measure(),
            vcf.ag_count(), vcf.other_mutations_count(), ans.size())

        if args.plot:
            name = os.path.basename(args.vcf_file).split('.')[0]
            bplt = BenchmarkPlot('plot_' + name)
            bplt.plot_accuracy(lab=str(vcf.vcf_name()),
                               recall=int(bench.recall()),
                               precision=int(bench.recall()))
            
if __name__ == '__main__':
    run()
    
