from __future__ import division
import vcf
import os.path
import re
import csv
import ConfigParser
import utils
from urllib2 import Request, urlopen, URLError
from collections import Counter

class DarnedDataGenerator(object):
    '''
    DarnedDataGenerator provides to prepare data that are used for benchmarking test.
    '''
    
    def __init__(self, species=None):
        __species = {
            'human_hg19':'http://darned.ucc.ie/static/downloads/hg19.txt',
            'human_hg18':'http://darned.ucc.ie/static/downloads/hg18.txt',
            'mice_mm9':'http://darned.ucc.ie/static/downloads/mm9.txt',
            'mice_mm10':'http://darned.ucc.ie/static/downloads/mm10.txt',
            'fly_dm3':'http://darned.ucc.ie/static/downloads/dm3.txt',
        }
        for k in __species:
            if k == species:
                self.species = species
                break;
        else: raise RuntimeError, "Given [%s] is not valid species name" % (species)

        if self.species is not None:
            self.filename = "".join([self.species, '.txt'])
            self.url = __species[self.species]
            self.saved_abs_path = utils.find_app_root() + '/data/'
            name, _ = os.path.splitext(self.filename)
            self.out_name = self.saved_abs_path + 'darned_' + name + '.csv'
            
    def fetch_darned(self):
        '''
        Fetch specify raw dataest from darned.ucc.ie/static/downloads/ into APP_ROOT/data,
        species name must be given, and acceptable type is defined as:
        human_hg18/hg19, mice_mm9/mm10, fly_dm3
        '''
        
        if os.path.isfile(self.saved_abs_path + self.filename):
            print "%s is already exist" % (self.filename)
            return False
        
        req = Request(self.url)
        try:
            response = urlopen(req)
        except URLError, e:
            if hasattr(e, 'reason'):
                print 'We failed to reach a server.'
                print 'Reason: ', e.reason
                return False
            elif hasattr(e, 'code'):
                print 'The server couldn\'t fulfill the request.'
                print 'Error code: ', e.code
                return False
        else:
            # works fine
            if not os.path.isdir(self.saved_path):
                os.makedirs(self.saved_path)
                print "Create directories [%s]" % (self.saved_path)
                
            print "Dowloading [%s] from [%s] ..." % (self.filename, self.url)
            with open(self.saved_path + self.filename, "w") as fout:
                fout.write(response.read())
            return True

    def darned_to_csv(self):
        '''
        Converting darned raw datafile to csv,
        the data that fetched from darned.ucc.ie/static/downloads/*.txt is given.
        >>> path_to_data = hg19.txt
        >>> darned_to_csv(path_to_data)
        Generate csv file into the APP_ROOT/data
        '''
        
        if not os.path.isfile(self.saved_path + self.filename):
            raise RuntimeError, 'Darned of %s is  not found' % (self.filename)

        if not os.path.isdir(self.saved_path):
            print "Create data dir"
            os.makedirs(self.saved_path)
        
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
    DarnedReader generates the subset of DARNED db
    >>> db = DarnedReader(sp='human_hg19', source='Brain', db='Path_to_Darned_DB')
    Returns list of subset of darned db
    '''
    
    def __init__(self, sp=None, source=None):
        if sp is None:
            raise RuntimeError, "Species name must be given"
        else: self.__sp = sp
        if source is None:
            self.__source = 'All'
        else:
            self.__source = source.upper()
            
        self.__darned_path = utils.find_app_root()+ '/data/'+ self.__sp+ '.csv'
        self.db = self.__generate_darned_set()
        
    def __str__(self):
        return "<%s.%s>" % (self.__class__.__name__)
        
    def __generate_darned_set(self):
        # Store selected records
        if self.__source is not None:
            selected = []

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
        elif self.__source is 'all' or self.__source is 'All':
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
        ''' given species name '''
        return self.__sp

    def path(self):
        ''' absolute path to Darned database file'''
        return os.path.abspath(self.__darned_path)
        
    def db_name(self):
        '''Darned db name'''
        return os.path.basename(self.__darned_path)
                 
    def size(self):
        '''number of the Darned entories'''
        return self.__db_size

        
class VCFReader(object):
    def __init__(self, filename):
        self.vcf = filename
        self.db = self.__generate_vcf_set()

    def __generate_vcf_set(self):
        ''' generate_vcf_set(self) -> list, returns the accumulated vcf'''
        vcf_recs = []
        vcf_reader = vcf.Reader(open(self.vcf, 'r'))
        self.count = 0
        self.substitutions = Counter()
        
        for rec in vcf_reader:
            types = str(rec.REF) + '-to-' + 'or'.join([str(i) for i in rec.ALT])
            self.substitutions[types] += 1
            mod_chr = re.sub(r'^chr', '', rec.CHROM, 1)
            vcf_recs.append(mod_chr+ ':'+ str(rec.POS))
            self.count += 1
        return vcf_recs
        
    def size(self):
        '''number of entory of the parsed vcf records'''
        return self.count

    def vcf_name(self):
        return os.path.basename(self.vcf)

    def editing_types(self):
        return self.substitutions
        
    def ag_count(self):
        return self.substitutions.get('A-to-G')
        
    def other_mutations_count(self):
        i = 0
        for k in self.substitutions:
            if not k == 'A-to-G':
                i += self.substitutions[k]
        return i
        
class Benchmark(object):
    '''
    >>> darned_db = DarnedReader(sp='human_hg19', source='Brain', db='Path_to_Darned_DB')
    >>> editing_db = VCFReader(filename)
    >>> bench = Benchmark(answer=darned_db, predict=candidate_db)
    '''
    
    def __init__(self, answer=None, predict=None):
        self.answer = set(answer)
        self.predict = set(predict)
        self.intersect = self.answer.intersection(self.predict)

    def __str__(self):
        return "DB[%d], Predict[%d]\n" % (len(self.answer), len(self.predict))

    def precision(self):
        precision = 0
        try:
            precision = len(self.intersect)/len(self.predict)
        except ZeroDivisionError:
            pass
        finally:
            return precision
            
    def recall(self):
        recall = 0
        try:
            recall = len(self.intersect)/len(self.answer)
        except ZeroDivisionError:
            pass
        finally:
            return recall
            
    def f_measure(self):
        precision = self.precision()
        recall = self.recall()
        
        f = 0
        try:
            f = (2*recall*precision)/(recall+precision)
        except ZeroDivisionError:
            pass
        finally:
            return f
