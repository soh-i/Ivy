from __future__ import division
import vcf
import os.path
import re
import csv
from Ivy.utils import Utils
from urllib2 import (
    Request,
    urlopen,
    URLError,
    )
from collections import Counter
import inspect

__program__ = 'benchmark'
__author__ = 'Soh Ishiguro <yukke@g-language.org>'
__license__ = ''
__status__ = 'development'


class DarnedDataGenerator(object):
    '''
    DarnedDataGenerator provides to prepare data that are used for benchmarking test.
    Args:
     species=(string)
    Returns:
     Darned db object
    Example:
     >>> ddg = DarnedDataGenerator(species=human_hg18)
    Raises:
     ValueError: when invalid species name was given
    '''
    
    def __init__(self, species=None):
        __species = {
            'human_hg19': 'http://darned.ucc.ie/static/downloads/hg19.txt',
            'human_hg18': 'http://darned.ucc.ie/static/downloads/hg18.txt',
            'mice_mm9': 'http://darned.ucc.ie/static/downloads/mm9.txt',
            'mice_mm10': 'http://darned.ucc.ie/static/downloads/mm10.txt',
            'fly_dm3': 'http://darned.ucc.ie/static/downloads/dm3.txt',
        }
        for k in __species:
            if k == species:
                self.species = species
                break;
        else:
            raise ValueError('Given species name: \'{sp}\' is invalid, {sps} are acceptable'.format(
                sp=species, sps="/".join([_ for _ in __species])))

        if self.species is not None:
            self.filename = ''.join([self.species, '.txt'])
            self.url = __species[self.species]
            self.saved_abs_path = Utils.find_app_root() + '/data/'
            name, _ = os.path.splitext(self.filename)
            self.out_name = self.saved_abs_path + 'darned_' + name + '.csv'
            
    def fetch_darned(self):
        '''
        Fetch specify raw dataest from darned.ucc.ie/static/downloads/ into the APP_ROOT/data.
        Create './data' directory if APP_ROOT/data is not found.
        Args:
         self
        Returns:
         bool
        Raises:
         URLError: could not connect Darned server
        '''
        
        if os.path.isfile(self.saved_abs_path + self.filename):
            sys.stderr.write("{f:s} is already exist".format(f=self.filename))
            return False
        
        req = Request(self.url)
        try:
            response = urlopen(req, timeout=10)
            
        except URLError, e:
            if hasattr(e, 'reason'):
                sys.stderr.write('We failed to reach a server due to {0}'.format(e.reason))
                raise URLError(", Could not connect " + req.get_full_url())
                
            elif hasattr(e, 'code'):
                sys.stderr.write('The server couldn\'t fulfill the request\n')
                sys.stderr.write('Error code: {e}'.format(e=e.code))
                raise URLError(", Could not connect " + req.get_full_url())
        else:
            # works fine
            if not os.path.isdir(self.saved_abs_path):
                os.makedirs(self.saved_abs_path)
                sys.stderr.write("Create directory into {path:s}\n".format(path=self.saved_abs_path))
                
            sys.stderr.write("Dowloading {filename:s} from {url:s} ...\n".format(
                filename=self.filename, url= self.url))
            
            with open(self.saved_abs_path + self.filename, "w") as fout:
                fout.write(response.read())
            return True

    def darned_to_csv(self):
        '''
        Converting darned raw datafile to csv, and generate csv file into the APP_ROOT/data.
        Args:
         self
        Returns:
         bool
        Raises:
         ValueError: when parsing error
        Exmples:
         >>> path_to_data = hg19.txt
         >>> darned_to_csv(path_to_data)
        '''
        
        if not os.path.isfile(self.saved_abs_path + self.filename):
            raise RuntimeError, 'Darned of %s is not found' % (self.filename)

        if not os.path.isdir(self.saved_abs_path):
            print "Create data dir"
            os.makedirs(self.saved_abs_path)
        
        if os.path.isfile(self.out_name):
            print "%s is already exisit" % (self.out_name)
            return False
        
        reader = csv.reader(open(self.saved_abs_path+self.filename, 'r'), delimiter="\t", quotechar="|")
        out = open(self.out_name, 'w')
        try:
            line_n = 0
            for row in reader:
                line_n += 1
                source = row[8]
                if len(source):
                    mod = source.replace(r';', ',').replace(r',', ';').\
                          replace(r'; ', ';').replace(r' ', '_').replace(r'_T', 'T')
                    out.write(",".join(row[:8]) + ",")
                    out.write(mod.upper() + ",")
                    out.write(",".join(row[9:]) + "\n")
        except:
            raise ValueError, 'Parsing error at line No.[%d]' % (line_n)
        finally:
            out.close()

class DarnedReader(object):
    '''
    DarnedReader generates subset of DARNED db.
    Args:
     sp=species name, source=speify tissue or cell line name
     Do not use source params, store to all records by the default settings.
    Example:
     >>> dr = DarnedReader(sp='human_hg19', source='Brain', db='Path_to_Darned_DB')
    Attributes:
     dr.db: array of subset of darned db
    '''
    
    def __init__(self, sp=None, source=None):
        if sp is None:
            raise RuntimeError, "Species name must be given"
        else:
            self.__sp = sp
            
        if source is None or len(source) == 0:
            self.__source = 'ALL'
        else:
            self.__source = source.upper()
            
        self.__darned_path = Utils.find_app_root()+ '/data/darned_'+ self.__sp+ '.csv'
        self.db = self.__generate_darned_set()
        
    def __str__(self):
        return "<%s.%s>" % (self.__class__.__name__)
        
    def __generate_darned_set(self):
        '''
        Generate darned db object stored as array
        '''
        # Store selected records
        if not self.__source == 'ALL':
            selected = []
            if os.path.isdir(self.__darned_path):
                raise IOError, "[%s] is directory, not csv file"
        
            with open(self.__darned_path, 'r') as f:
                for line in f:
                    if not line.startswith('chrom'):
                        _data = line.split(",")
                        _chrom = _data[0]
                        _pos = _data[1]
                        _darned_source = _data[8]
                        _each_source = _darned_source.split(";")
                        if any([_ for _ in _each_source if _ == self.__source]):
                            selected.append(':'.join([_chrom, _pos, self.__source]))
                            
                self.__db_size = len(selected)
                return selected
                
        # Store all Darned records (default)
        elif self.__source == 'ALL':
            darned_list = []
            with open(self.__darned_path, 'r') as f:
                for line in f:
                    if not line.startswith('chrom'):
                        _data = line.split(',')
                        _chrom = _data[0]
                        _pos = _data[1]
                        _darned_source = _data[8]
                        darned_list.append(':'.join([_chrom, _pos, _darned_source]))
                self.__db_size = len(darned_list)
                return darned_list
                
        elif not self.__sp:
            raise RuntimeError, 'Given species name[%s] is not valid' % (self.__sp)

    def sp(self):
        '''
        Returns:
         tuple: species and genome version
        '''
        (sp, ver) = self.__sp.split("_")
        return sp, ver

    def path(self):
        '''
        Returns:
         string: absolute path to Darned database file
        '''
        return os.path.abspath(self.__darned_path)
        
    def db_name(self):
        '''
        Returns:
         string: Darned db name
        '''
        return os.path.basename(self.__darned_path)
                 
    def size(self):
        '''
        Returns:
         int: number of the Darned entories
        '''
        return self.__db_size

        
class VCFReader(object):
    '''
    VCFReader class provides that list of VCF file within utils methods
    Args:
     filename(string): filename of VCF
    Returns:
     list: vcf files
    Attributes:
     db(list): Stored VCF list
    Examples:
     >>> vr = VCFReader(path_to_vcf_file)
    '''
    
    def __init__(self, filename):
        self.__vcf = filename
        self.db = self.__generate_vcf_set()

    def __generate_vcf_set(self):
        _vcf_recs = []
        
        if os.path.isdir(self.__vcf):
            raise IOError, "[%s] is directory, not csv file"
        
        _vcf_reader = vcf.Reader(open(self.__vcf, 'r'))
        self.__substitutions = Counter()
        
        for rec in _vcf_reader:
            _types = str(rec.REF) + '-to-' + 'or'.join([str(_) for _ in rec.ALT])
            self.__substitutions[_types] += 1
            _mod_chr = re.sub(r'^chr', '', rec.CHROM, 1)
            _vcf_recs.append(_mod_chr+ ':'+ str(rec.POS))
        self.__size = len(_vcf_recs)
        return _vcf_recs
        
    def size(self):
        '''
        Returns:
         int: number of entory of the parsed vcf records
        '''
        return self.__size

    def vcf_name(self):
        '''
        Returns:
         string: filename of vcf
        '''
        return os.path.basename(self.__vcf)

    def editing_types(self):
        '''
        Returns:
         collection object: All type of the base substitutions
        '''
        return self.__substitutions
        
    def ag_count(self):
        '''
        Returns:
         int: A-to-G editing alone
        '''
        return self.__substitutions.get('A-to-G')

    def target_count(self, types):
        '''
        Args:
         types(string): specify type of base substitution
        Returns:
         int: base substitution count
        '''
        return self.__substitutions.get(types)
        
    def other_mutations_count(self):
        '''
        Returns:
         collection object: substitution count except A-to-G editing
        '''
        i = 0
        for k in self.__substitutions:
            if not k == 'A-to-G':
                i += self.__substitutions[k]
        return i
        
    
class __CSVReader(object):
    '''
    CSVReader class provides to generate array of CSV file
    Examples:
     >>> csv = VCFReader(path_to_csv_file)
     >>> csv.db
    '''
    
    def __init__(self, filename):
        self.__filename = filename
        self.db = self.__generate_csv_set()
        
    def __generate_csv_set(self):
        csv_recs = []
        if os.path.isdir(self.__filename):
            raise IOError, "[%s] is directory, not csv file"
            
        __line__ = 1
        with open(self.__filename) as f:
            for line in f:
                #if not line.startswith("track") \
                #  and not line.startswith('#') \
                #  and not line.startswith("Chromosome") \
                #  and not line.startswith("Ch,") \
                #  and not line.startswith("Arm,"):
                if __line__ != 1 \
                   and not line.startswith('#') \
                   and not line.startswith('Chr,'):
                    rec = line.split(',')
                    if rec[0].startswith('chr'):
                        _chrom = re.sub(r'^chr', '', rec[0], 1)
                        str(_chrom)
                    elif len(rec[0]):
                        _chrom = str(rec[0])
                    if rec[1].find(','):
                        _pos = rec[1].replace(',', '')
                        str(_pos)
                    csv_recs.append(_chrom + ':' + _pos)
                __line__ += 1
                
        self.__size = len(csv_recs)
        return csv_recs

    def size(self):
        return self.__size

    def name(self):
        return os.path.basename(self.__filename)

        
class Benchmark(object):
    '''
    Benchmark calculates precision, recall and F-measure from given data set,
    those metrics are defined as follows:
    precision = TP/(TP+FP), recall = TP/(TP+FN),
    F-measure = 2*Precision*Recall/(Precision+Recall)
    Args:
     answer(set), predict(set)
    Examples:
     >>> >>> bench = Benchmark(answer=darned_db, predict=candidate_db)
    Attributes:
     answer(set): set of answer sites from Darned
     predict(set): set of predicted sites
     intersect(set): intersection between anseer and predict
    '''
    
    def __init__(self, answer=None, predict=None):
        if not isinstance(answer, list):
            raise TypeError, "[%s] is given, data must be list alone" % (type(answer))
        elif not isinstance(predict, list):
            raise TypeError, "[%s] is given, data must be list alone" % (type(predict))
            
        # remove string(tissue/sample info) except chromosome and position
        self.answer = set([":".join(_.split(":")[:2]) for _ in answer])
        self.predict = set([":".join(_.split(":")[:2]) for _ in predict])
        
        if len(self.answer) == 0:
            raise ValueError, 'Answer data set has no entory'
        elif len(self.predict) == 0:
            raise ValueError, 'Candidate data set has no entory'
                
        self.intersect = self.answer.intersection(self.predict)

    def __str__(self):
        return "Answer set[%d], Candidate set[%d]\n" % (len(self.answer), len(self.predict))

    def precision(self):
        '''
        Returns:
         float: precision
        '''
        try:
            _precision = len(self.intersect)/len(self.predict)
            return _precision
        except ZeroDivisionError:
            _precision = 0
        finally:
            return _precision
            
    def recall(self):
        '''
        Returns:
         float: recall
        '''
        try:
            _recall = len(self.intersect)/len(self.answer)
            return _recall
        except ZeroDivisionError:
            _recall = 0
        finally:
            return _recall
            
    def f_measure(self):
        '''
        Returns:
         float: F-measure
        '''
        _precision = self.precision()
        _recall = self.recall()
        try:
            _f = 2*_recall*_precision/(_recall+_precision)
            return _f
        except ZeroDivisionError:
            _f = 0
        finally:
            return _f
