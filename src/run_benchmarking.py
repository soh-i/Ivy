import os.path
import sys
from Ivy.commandline.parse_benchmarking_opts import parse_bench_opts
from Ivy.benchmark.plot import BenchmarkPlot
from Ivy.benchmark.benchmark import *
import Ivy.utils

__program__ = 'ivy_benchmark'
__author__ = 'Soh Ishiguro <yukke@g-language.org>'
__license__ = ''
__status__ = 'development'

def run():
    args = parse_bench_opts()
    _prepare_required_data(args.sp)

    if args.vcf_file:
        result = _call_bench(args.vcf_file, sp=args.sp, source=args.source, mode='vcf')
    elif args.csv_file:
        result = _call_bench(args.csv_file, sp=args.sp, source=args.source, mode='csv')

    if args.out:
        _write_result(filename=args.out, content=result, is_file=True)
    elif not args.out:
        _write_result(content=result, is_file=False)
        
    if args.plot:
        _plot(precisions, recalls, args.csv_file)
    
def _prepare_required_data(species):
    try:
        gen = DarnedDataGenerator(species)
    except DarnedDataGeneratorValueError as e:
        raise SystemExit('[{cls}]: given \'{sp}\' is not valid name, {sps} is only valid name'.format(
            cls=e.__class__.__name__, sp=e.sp, sps=e.sps))
        
    darned_raw_file = gen.saved_abs_path + gen.filename
    if not os.path.isfile(darned_raw_file):
        sys.stderr.write("No {f:s} file\n".format(f=darned_raw_file))
        sys.stderr.write("Fetching {sp:s} from Darned DB...\n".format(sp=species))
        gen.fetch_darned()
        
    darned_parsed_csv = gen.out_name
    if not os.path.isfile(darned_parsed_csv):
        sys.stderr.write("Parsing Darned db...\n")
        try:
            gen.darned_to_csv()
        except DarnedDataGeneratorParseError as e:
            raise SystemExit('[{cls}]: {e}'.format(cls=e.__class__.__name__, e=e))

def _call_bench(vcf_files, sp=None, source=None, mode=None):
    precisions, recalls, f_measures = [0, 0, 0]
    content = str()
    ans = DarnedReader(sp=sp, source=source)
    
    for v in vcf_files:
        if mode == 'vcf':
            pred = VCFReader(v)
        elif mode == 'csv':
            pred = __CSVReader()
        try:
            bench = Benchmark(answer=ans.db, predict=pred.db)
        except BenchmarkIOException as e:
            raise SystemExit('[{0}]: {1}'.format(e.__class__.__name__, e))
            
        p = bench.precision()
        r = bench.recall()
        f = bench.f_measure()
        content += (
                '{species:s},{source:s},{db_name:s},{vcf_file:s},{precision:f},{recall:f},{f_measure:f},{ag_count:d},{other_count:d},{ans_size:d}\n'
                .format(
                    species=ans.sp()[0],
                    source=source,
                    db_name=ans.db_name(),
                    vcf_file=pred.name(),
                    precision=p,
                    recall=r,
                    f_measure=f,
                    ag_count=pred.ag_count(),
                    other_count=pred.other_mutations_count(),
                    ans_size=ans.size()))
    return content
        
def _write_result(is_file=False, **data):
    header = "Species,Source,DB,VCF,Precision,Recall,F-measure,AGs,Others,AnsCount\n"
    
    if is_file is False:
        sys.stdout.write(header)
        sys.stdout.write(data['content'])
    elif is_file is True:
        f = open(data['filename'], 'w')
        f.write(header)
        f.write(data['content'])
        f.close()
    
def _use_csv():
    print "Species,Source,DB,CSV,Precision,Recall,F-measure,PredCount,AnsCount"

    try:
        ans = DarnedReader(sp=args.sp, source=args.source)
    except TypeError as e:
        raise SystemExit('[{cls}]: {err}'.format(cls=e.__class__.__name__, err=e))
    except BenchmarkIOException as e:
        raise SystemExit('[{cls}]: {err}'.format(cls=e.__class__.__name__, err=e))
            
    precisions = []
    recalls = []
    f_measures = []
        
    for c in args.csv_file:
        csv = __CSVReader(c)
        bench = Benchmark(answer=ans.db, predict=csv.db)
        p = bench.precision()
        r = bench.recall()
        f = bench.f_measure()
            
        print "%s,%s,%s,%s,%f,%f,%f,%d,%d" % (
            ans.sp()[0],
            args.source,
            ans.db_name(),
            csv.name(),
            p, r, f,
            csv.size(),
            ans.size())

        precisions.append(float(p))
        recalls.append(float(r))
        f_measures.append(float(f))

def __plot(p, r, labs):
    if isinstance(p, list) and isinstance(r, list) and isinstance(labs, list):
        names = [os.path.basename(_).split('.')[0] for _ in labs]
        #bplt = BenchmarkPlot('plot_' + ','.join(names), "human")
        bplt = BenchmarkPlot('plot_' + "test", "human")
        bplt.plot_accuracy(lab=names, recall=r, precision=p)
        return True
    else:
        return False
