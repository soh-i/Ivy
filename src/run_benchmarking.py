from Ivy.benchmark.benchmark import (
    DarnedDataGenerator,
    DarnedReader,
    VCFReader,
    Benchmark,
    __CSVReader,
    )
from Ivy.benchmark.plot import BenchmarkPlot
from Ivy.version import __version__
import Ivy.utils
import argparse
import os.path

__program__ = 'ivy_benchmark'
__author__ = 'Soh Ishiguro <yukke@g-language.org>'
__license__ = ''
__status__ = 'development'

def run():
    desc = "Benchmarking test for detected RNA editing sites based on HTSeq data."
    
    parser = argparse.ArgumentParser(description=desc)
    group = parser.add_mutually_exclusive_group(required=True)
    
    group.add_argument('--vcf',
                       dest='vcf_file',
                       action='store',
                       nargs='+',
                       help='set VCF file [required]'
                   )
    group.add_argument('--csv',
                       dest='csv_file',
                       action='store',
                       nargs='+',
                       help='set CSV file [required]'
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
                        help='set species and genome version (eg. human_hg19) [required]'
                    )
    parser.add_argument('--plot',
                       required=False,
                       default=False,
                       action='store_true',
                       help='plot benchmarking stats (default:False)',
                   )
    parser.add_argument('--version', action='version', version=__version__)
    
    args = parser.parse_args()
    gen = DarnedDataGenerator(species=args.sp)
    
    darned_raw_file = gen.saved_abs_path + gen.filename
    if not os.path.isfile(darned_raw_file):
        print "Fetching from darned..."
        gen.fetch_darned()
        
    darned_parsed_csv = gen.out_name
    if not os.path.isfile(darned_parsed_csv):
        print "parsing darned..."
        gen.darned_to_csv()

    # use VCF files
    if args.vcf_file and args.sp:
        print "Species,DB,VCF,Precision,Recall,F-measure,AGs,Others,AnsCount"
        
        ans = DarnedReader(sp=args.sp, source=args.source)
        precision = []
        recall = []
        f_measure = []
        for v in args.vcf_file:
            vcf = VCFReader(v)
            bench = Benchmark(answer=ans.db, predict=vcf.db)
            p = bench.precision()
            r = bench.recall()
            f = bench.f_measure()
            
            print "%s,%s,%s,%f,%f,%f,%d,%d,%d" % (
                ans.sp()[0],
                ans.db_name(),
                vcf.vcf_name(),
                p, r, f,
                vcf.ag_count(),
                vcf.other_mutations_count(),
                ans.size())
            
            precision.append(float(p))
            recall.append(float(r))
            f_measure.append(float(f))
            
        if args.plot:
            names = [os.path.basename(_).split('.')[0] for _ in args.vcf_file]
            bplt = BenchmarkPlot('plot_' + ','.join(names))
            bplt.plot_accuracy(lab=args.vcf_file, recall=r, precision=p)
            
    # use CSV files
    elif args.csv_file and args.sp:
        print "Species,DB,CSV,Precision,Recall,F-measure,AnsCount"

        ans = DarnedReader(sp=args.sp, source=args.source)
        precision = []
        recall = []
        f_measure = []
        
        for c in args.csv_file:
            csv = __CSVReader(c)
            bench = Benchmark(answer=ans.db, predict=csv.db)
            p = bench.precision()
            r = bench.recall()
            f = bench.f_measure()

            print "%s,%s,%s,%f,%f,%f,%d" % (
                ans.sp()[0],
                ans.db_name(),
                csv.name(),
                p, r, f,
                ans.size())

            precision.append(float(p))
            recall.append(float(r))
            f_measure.append(float(f))
        
        if args.plot:
            names = [os.path.basename(_).split('.')[0] for _ in args.csv_file]
            bplt = BenchmarkPlot('plot_' + ','.join(names))
            bplt.plot_accuracy(lab=args.csv_file, recall=r, precision=p)
            
