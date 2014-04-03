import os.path
import sys
import datetime
from Ivy.cli.edit_bench_cli_opts import parse_bench_opts
from Ivy.analysis_settings import EDIT_BENCH_SETTINGS
from Ivy.edit_bench.plot import BenchmarkPlot
from Ivy.edit_bench.benchmark import (
    DarnedDataGenerator,
    DarnedReader,
    VCFReader,
    Benchmark,
    CSVReader,
)

__program__ = 'edit_bench'
__author__ = 'Soh Ishiguro <yukke@g-language.org>'
__license__ = 'GPL v2'
__status__ = 'development'


class AppEditBench(object):
    def __init__(self):
        self.args = parse_bench_opts()
        self.__prepare_required_data(self.args.sp)

    def run(self):
        # Input format as VCF
        if self.args.vcf_file and self.args.sp:
            if self.args.plot:
                result = self._call_bench(self.args.vcf_file, sp=self.args.sp, source=self.args.source,
                                          mode='vcf', plt=True, labs=self.args.vcf_file)
            else:
                result = self._call_bench(self.args.vcf_file, sp=self.args.sp, source=self.args.source,
                                          mode='vcf')
        # Input format as CSV
        elif self.args.csv_file and self.args.sp:
            if self.args.plot:
                result = self._call_bench(self.args.csv_file, sp=self.args.sp, source=self.args.source,
                                     mode='csv', plt=True, labs=self.args.csv_file)
            else:
                # Make p-r plot
                result = self._call_bench(self.args.csv_file, sp=self.args.sp, source=self.args.source,
                                     mode='csv')
        # Output filename
        if self.args.out:
            self._write_result(filename=self.args.out, content=result, is_file=True)
        elif not self.args.out:
            self._write_result(content=result, is_file=False)
    
    def __prepare_required_data(self, species):
        try:
            gen = DarnedDataGenerator(species)
        except ValueError as e:
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

    def _call_bench(self, files, sp=None, source=None, mode=None, plt=False, labs=[]):
        precisions, recalls, f_measures = [], [], []
        content = str()
        ans = DarnedReader(sp=sp, source=source)
    
        for f in files:
            if mode == 'vcf':
                pred = VCFReader(f)
            elif mode == 'csv':
                pred = CSVReader(f)
            else:
                raise ValueError("Valid input filename is ONLY csv or vcf")
            
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
            
        # plot data
        if len(precisions) > 0 and len(recalls) > 0:
            try:
                self.__plot(precisions, recalls, labs)
            except TypeError as e:
                raise SystemExit(e)
                
        return content
        
    def _write_result(self, is_file=False, **data):
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
            with open(data['filename'], 'w') as f:
                f.write(header)
                f.write(data['content'])
            
    def __plot(self, precision, recall, labs):
        if isinstance(precision, list) and isinstance(recall, list) and isinstance(labs, list):
            names = [os.path.basename(_).split('.')[0] for _ in labs]
            d = datetime.datetime.today()
            time = str(d.now()).split(" ")[1]
            bplt = BenchmarkPlot('plot_' + time,  self.args.sp)
            bplt.plot_accuracy(lab=names, recall=recall, precision=precision)
            return True
        else:
            raise TypeError("[Error] Input data type must be \'list\' to plot data")

            
if __name__ == '__main__':
    app = AppEditBench()
    app.run()
    
