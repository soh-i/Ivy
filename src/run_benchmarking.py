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

    if args.vcf_file and args.sp:
        result = _call_bench(args.vcf_file, sp=args.sp, source=args.source)
        
    header = "Species,Source,DB,VCF,Precision,Recall,F-measure,AGs,Others,AnsCount\n"
    if args.out:
        _write_result(args.out, header, result)
    elif not args.out:
        sys.stdout.write(header)
        sys.stdout.write(result)
        
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

def _call_bench(vcf_files, sp=None, source=None):
    precision, recall, f_measure = [0, 0, 0]
    content = str()
    ans = DarnedReader(sp=sp, source=source)
    for v in vcf_files:
        vcf = VCFReader(v)
        try:
            bench = Benchmark(answer=ans.db, predict=vcf.db)
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
                    vcf_file=vcf.vcf_name(),
                    precision=p,
                    recall=r,
                    f_measure=f,
                    ag_count=vcf.ag_count(),
                    other_count=vcf.other_mutations_count(),
                    ans_size=ans.size()))
        return content
        
def _write_result(filename, header, content):
    f = open(filename, 'w')
    f.write(header)
    f.write(content)
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
