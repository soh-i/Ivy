#!/usr/bin/env python

from Ivy.benchmark import DarnedDataGenerator, DarnedReader, VCFReader, Benchmark
import Ivy.utils
import argparse
import os.path

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
                        help='use specific sample/tissue/cell line, default:All'
                    )
    parser.add_argument('--sp',
                        required=True,
                        dest='sp',
                        action='store',
                        help='set species name and genome version, ex: human_hg19 [required]'
                    )
    parser.add_argument('--version', action='version', version=version)
    
    args = parser.parse_args()
    generator = DarnedDataGenerator(species=args.sp)

    darned_raw_file = generator.saved_path + generator.filename
    if not os.path.isfile(darned_raw_file):
        generator.fetch_darned()

    darned_parsed_csv = generator.out_name
    if not os.path.isfile(darned_parsed_csv):
        generator.darned_to_csv()

    if args.vcf_file and args.sp and args.source:
        ans = DarnedReader(sp=args.sp, source=args.source)
        vcf = VCFReader(args.vcf_file)
        bench = Benchmark(answer=ans.db, predict=vcf.db)
        
        print "Species:%s\tDB:%s\tVCF:%s\tPrecision:%f\tRecall:%f\tF-measure:%f\tAGs:%d\tOthers:%d\tAnsCount:%d\n" % (
            ans.sp, ans.db_name(), vcf.vcf_name(),
            bench.precision(), bench.recall(), bench.f_measure(),
            vcf.ag_count(), vcf.other_mutations_count(), vcf.size()),
        
