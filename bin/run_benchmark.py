#!/usr/bin/env python

from Ivy.benchmark import DarnedReader, VCFReader, Benchmark
import argparse

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
    parser.add_argument('--sp',
                        required=True,
                        dest='species',
                        action='store',
                        help='set species name and genome version, ex: human_hg19 [required]'
                    )
    parser.add_argument('--db',
                        required=False,
                        dest='db_file',
                        action='store',
                        help='set db name [not implemented]'
                    )
    parser.add_argument('--config',
                        required=False,
                        dest='conf_file',
                        action='store',
                        help='set config file [not implemented]'
                    )
    parser.add_argument('--version', action='version', version=version)
    args = parser.parse_args()
    
    if args.vcf_file and args.species:
        ans = DarnedReader(sp=args.species)
        vcf = VCFReader(args.vcf_file)
        bench = Benchmark(answer=ans.db, predict=vcf.db)
        
        print "Species:%s\tDB:%s\tVCF:%s\tPrecision:%f\tRecall:%f\tF-measure:%f\tAGs:%d\tOthers:%d\tAnsCount:%d\n" % (
            ans.sp(), ans.db_name(), vcf.vcf_name(),
            bench.precision(), bench.recall(), bench.f_measure(),
            vcf.ag_count(), vcf.other_mutations_count(), ans.count()),
        
