import os.path
import sys
from Ivy.commandline.parse_benchmarking_opts import parse_bench_opts
from Ivy.benchmark.plot import BenchmarkPlot
from Ivy.benchmark.benchmark import (
    DarnedDataGenerator,
    DarnedDataGeneratorValueError,
    DarnedDataGeneratorParseError,
    DarnedReader,
    VCFReader,
    Benchmark,
    BenchmarkIOException,
    __CSVReader,
)
import Ivy.utils

__program__ = 'ivy_benchmark'
__author__ = 'Soh Ishiguro <yukke@g-language.org>'
__license__ = ''
__status__ = 'development'

def run():
    '''
    Run benchmarking test when call this function
    '''
    args = parse_bench_opts()
    _prepare_required_data(args.sp)

    # Input format is vcf
    if args.vcf_file and args.sp:
        if args.plot:
            result = _call_bench(args.vcf_file, sp=args.sp, source=args.source,
                                 mode='vcf', plt=True, labs=args.vcf_file)
        else:
            result = _call_bench(args.vcf_file, sp=args.sp, source=args.source, mode='vcf')
    # Input format is csv
    elif args.csv_file and args.sp:
        if args.plot:
            result = _call_bench(args.csv_file, sp=args.sp, source=args.source,
                                 mode='csv', plt=True, labs=args.csv_file)
        else:
            result = _call_bench(args.csv_file, sp=args.sp, source=args.source, mode='csv')
        
    # Output
    if args.out:
        _write_result(filename=args.out, content=result, is_file=True)
    elif not args.out:
        _write_result(content=result, is_file=False)
    
def _prepare_required_data(species):
    '''
    prepares and generates the required data files to test benchmaring
    
    Args:
     species(string): must be given species name
    
    Raises:
     DarnedDataGeneratorValueError: Invalid species name was given
     DarnedDataGeneratorParseError: Faild to parse csv from Darned row data file
    '''
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

def _call_bench(files, sp=None, source=None, mode=None, plt=False, labs=[]):
    '''
    Wrapper of benchmark class

    Args:
     files(string): vcf/csv file name
     sp(string): species name
     source(string): exploring specified sample/tissues
     mode(string): select input file fomrat
     plt(dict): plot data or not

    Raises:
     ValueError: invalid file format as mode args
     SystemExit: fetal error when create benchmarking instance
    
    Returns:
     content(string): all of the benchmarking results
    '''
    precisions, recalls, f_measures = [], [], []
    content = str()
    ans = DarnedReader(sp=sp, source=source)
    
    for f in files:
        if mode == 'vcf':
            pred = VCFReader(f)
        elif mode == 'csv':
            pred = __CSVReader(f)
        else:
            raise ValueError("Valid input filename is csv or vcf alone")
            
        try:
            bench = Benchmark(answer=ans.db, predict=pred.db)
        except BenchmarkIOException as e:
            raise SystemExit('[{0}]: {1}'.format(e.__class__.__name__, e))
            
        p = bench.precision()
        r = bench.recall()
        f = bench.f_measure()
        content += (
                '{sp:s},{src:s},{db_name:s},{vcf_f:s},{precision:f},{recall:f},{f_measure:f},{p_size},{a_size:d}\n'
                .format(
                    sp=ans.sp()[0],
                    src=source,
                    db_name=ans.db_name(),
                    vcf_f=pred.name(),
                    precision=p,
                    recall=r,
                    f_measure=f,
                    p_size=pred.size(),
                    a_size=ans.size()))
        
        if plt is True and len(labs) > 0:
            precisions.append(p)
            recalls.append(r)
    if len(precisions) > 0 and len(recalls) > 0:
        __plot(precisions, recalls, labs)
        
    return content
        
def _write_result(is_file=False, **data):
    '''
    Write benchmarking result into the file or console
    
    Args:
     is_file(bool): print file or stdout
     data(dict): character of the benchmarking and header
    '''
    header = "Species,Source,DB,VCF,Precision,Recall,F-measure,AGs,Others,AnsCount\n"
    
    if is_file is False:
        sys.stdout.write(header)
        sys.stdout.write(data['content'])
    elif is_file is True:
        f = open(data['filename'], 'w')
        f.write(header)
        f.write(data['content'])
        f.close()

def __plot(p, r, labs):
    '''
    Plot benchmarking result
    
    Args:
     p(float): precison
     r(float): recall
     labs(list): labels of plot
    '''
    if isinstance(p, list) and isinstance(r, list) and isinstance(labs, list):
        names = [os.path.basename(_).split('.')[0] for _ in labs]
        #bplt = BenchmarkPlot('plot_' + ','.join(names), "human")
        bplt = BenchmarkPlot('plot_' + "test", "human")
        bplt.plot_accuracy(lab=names, recall=r, precision=p)
        return True
    else:
        return False
