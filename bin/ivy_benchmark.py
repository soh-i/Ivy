#!/usr/bin/env python

from Ivy.benchmark.benchmark import *
from Ivy.benchmark.plot import BenchmarkPlot
import Ivy.utils
import argparse
import os.path

__program__ = 'benchmark'
__author__ = 'Soh Ishiguro <yukke@g-language.org>'
__status__ = 'development'

if __name__ == '__main__':
    version = '0.0.1'
    desc = "Benchmarking test for detected RNA editing sites based on HTSeq data."
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument('--vcf',
                        required=True,
                        dest='vcf_file',
                        action='store',
                        help='set VCF file [required]'
                    )
    parser.add_argument('--source',
                        required=False,
                        dest='source',
                        action='store',
                        help='use specific sample/tissue/cell line (default:All)'
                    )
    parser.add_argument('--sp',
                        required=True,
                        dest='sp',
                        action='store',
                        help='set species and genome version (eg. ehuman_hg19) [required]'
                    )
    parser.add_argument('--plot',
                       required=False,
                       default=False,
                       action='store_true',
                       help='plot benchmarking stats (default:False)',
                   )
    parser.add_argument('--version', action='version', version=version)
    
    args = parser.parse_args()
    generator = DarnedDataGenerator(species=args.sp)
    
    darned_raw_file = generator.saved_abs_path + generator.filename
    if not os.path.isfile(darned_raw_file):
        generator.fetch_darned()
        
    darned_parsed_csv = generator.out_name
    if not os.path.isfile(darned_parsed_csv):
        generator.darned_to_csv()

        
    if args.vcf_file and args.sp and args.source:
        ans = DarnedReader(sp=args.sp, source=args.source)
        vcf = VCFReader(args.vcf_file)
        bench = Benchmark(answer=ans.db, predict=vcf.db)
        
        print "Species:%s,DB:%s,VCF:%s,Precision:%f,Recall:%f,F-measure:%f,AGs:%d,Others:%d,AnsCount:%d" % (
            ans.sp()[0], ans.db_name(), vcf.vcf_name(),
            bench.precision(), bench.recall(), bench.f_measure(),
            vcf.ag_count(), vcf.other_mutations_count(), ans.size())

        if args.plot:
            bplt = BenchmarkPlot("testplot")
            bplt.plot_accuracy(lab=str(vcf.vcf_name()),
                               recall=int(bench.recall()),
                               precision=int(bench.recall()))
            
    
        
        
