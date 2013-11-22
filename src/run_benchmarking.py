from argparse import ArgumentParser
import os.path
import sys
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

__program__ = 'ivy_benchmark'
__author__ = 'Soh Ishiguro <yukke@g-language.org>'
__license__ = ''
__status__ = 'development'

def run():
    desc = "Benchmarking test for detected RNA editing sites based on HTSeq data to evaluate detection params."
    
    parser = ArgumentParser(description=desc,
                            prog=__program__,
                            )
    group = parser.add_mutually_exclusive_group(required=True)
    
    group.add_argument('--vcf',
                       dest='vcf_file',
                       action='store',
                       nargs='+',
                       metavar='',
                       help='VCF file(s)',
                       )
    group.add_argument('--csv',
                       dest='csv_file',
                       action='store',
                       nargs='+',
                       metavar='',
                       help='CSV file(s), For debug mode',
                       )
    parser.add_argument('--source',
                        required=False,
                        default='All',
                        dest='source',
                        action='store',
                        metavar='',
                        help='specific sample/tissue/cell line [default: All]',
                        )
    parser.add_argument('--sp',
                        required=True,
                        dest='sp',
                        metavar='species',
                        action='store',
                        help='species and genome version (eg. human_hg19)',
                        )
    parser.add_argument('--plot',
                        required=False,
                        default=False,
                        action='store_true',
                        help='plot benchmarking stats [default: Off]',
                        )
    parser.add_argument('--out',
                        required=False,
                        action='store_true',
                        help='output name',
                        )
    parser.add_argument('--version',
                        action='version',
                        version=__version__
                        )

    args = parser.parse_args()
    gen = DarnedDataGenerator(species=args.sp)
    
    # check darned raw file
    darned_raw_file = gen.saved_abs_path + gen.filename
    if not os.path.isfile(darned_raw_file):
        sys.stderr.write("No {f:s} file\n".format(f=darned_raw_file))
        sys.stderr.write("Fetching {sp:s} from Darned DB...\n".format(sp=args.sp))
        gen.fetch_darned()
        
    darned_parsed_csv = gen.out_name
    if not os.path.isfile(darned_parsed_csv):
        sys.stderr.write("Parsing Darned db...\n")
        gen.darned_to_csv()

    # use VCF files
    if args.vcf_file and args.sp:
        print "Species,Source,DB,VCF,Precision,Recall,F-measure,AGs,Others,AnsCount"
        
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

            print '{species:s}\t{source:s}\t{db_name:s}\t{vcf_file:s}\t{precision:f}\t \
            {recall:f}\t{f_measure:f}\t{ag_count:d}\t{other_count:d}\t{ans_size:d}'.format(
                species=ans.sp()[0],
                source=args.source,
                db_name=ans.db_name(),
                vcf_file=vcf.vcf_name(),
                precision=p,
                recall=r,
                f_measure=f,
                ag_count=vcf.ag_count(),
                other_count=vcf.other_mutations_count(),
                ans_size=ans.size())
            
        if args.plot:
            names = [os.path.basename(_).split('.')[0] for _ in args.vcf_file]
            bplt = BenchmarkPlot('plot_' + ','.join(names))
            bplt.plot_accuracy(lab=args.vcf_file, recall=r, precision=p)
            
    # use CSV files, this is debug mode
    elif args.csv_file and args.sp:
        print "Species,Source,DB,CSV,Precision,Recall,F-measure,AnsCount"

        ans = DarnedReader(sp=args.sp, source=args.source)
        precisions = []
        recalls = []
        f_measures = []
        
        for c in args.csv_file:
            csv = __CSVReader(c)
            bench = Benchmark(answer=ans.db, predict=csv.db)
            p = bench.precision()
            r = bench.recall()
            f = bench.f_measure()

            print "%s,%s,%s,%s,%f,%f,%f,%d" % (
                ans.sp()[0],
                args.source,
                ans.db_name(),
                csv.name(),
                p, r, f,
                ans.size())

            precisions.append(float(p))
            recalls.append(float(r))
            f_measures.append(float(f))
        
        if args.plot:
            plot(precisions, recalls, args.csv_file)

def plot(p, r, labs):
    if isinstance(p, list) and isinstance(r, list) and isinstance(labs, list):
        names = [os.path.basename(_).split('.')[0] for _ in labs]
        #bplt = BenchmarkPlot('plot_' + ','.join(names), "human")
        bplt = BenchmarkPlot('plot_' + "test", "human")
        bplt.plot_accuracy(lab=names, recall=r, precision=p)
        return True
    else:
        return False
