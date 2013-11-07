from Ivy.benchmark.benchmark import (
    DarnedDataGenerator,
    DarnedReader,
    Benchmark,
    CSVReader,
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
    version = __version__
    desc = "For CSV data"
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument('--csv',
                        required=True,
                        dest='csv_file',
                        action='store',
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
    parser.add_argument('--version', action='version', version=version)
    
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
        
    if args.csv_file and args.sp:
        ans = DarnedReader(sp=args.sp, source=args.source)
        csv = CSVReader(args.csv_file)
        bench = Benchmark(answer=ans.db, predict=csv.db)
        
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
    
