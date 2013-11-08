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

__program__ = 'benchmark'
__author__ = 'Soh Ishiguro <yukke@g-language.org>'
__license__ = ''
__status__ = 'development'

class DarnedDataGenerator(object):
    '''
    DarnedDataGenerator provides to prepare data that are used for benchmarking test.
    >>> ddg = DarnedDataGenerator(species=human_hg18)
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
            self.saved_abs_path = Utils.find_app_root() + '/data/'
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
            response = urlopen(req, timeout=10)
        except URLError, e:
            if hasattr(e, 'reason'):
                print 'We failed to reach a server.'
                print 'Reason: ', e.reason
                raise URLError(", Could not connect " + req.get_full_url())
            elif hasattr(e, 'code'):
                print 'The server couldn\'t fulfill the request.'
                print 'Error code: ', e.code
                raise URLError(", Could not connect " + req.get_full_url())
        else:
            # works fine
            if not os.path.isdir(self.saved_abs_path):
                os.makedirs(self.saved_abs_path)
                print "Create directory into [%s]" % (self.saved_abs_path)
                
            print "Dowloading [%s] from [%s] ..." % (self.filename, self.url)
            with open(self.saved_abs_path + self.filename, "w") as fout:
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
    DarnedReader class generates subset of DARNED db.
    >>> dr = DarnedReader(sp='human_hg19', source='Brain', db='Path_to_Darned_DB')
    >>> dr.db
    Returns array of subset of darned db.
    Do not use source option, store to all records by the default settings.
    Acceptable type of sp argument is defined as human_hg18/hg19, mice_mm9/mm10, fly_dm3.
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
        # Store selected records
        if not self.__source == 'ALL':
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
        '''return tuple of species and genome version'''
        (sp, ver) = self.__sp.split("_")
        return sp, ver

    def path(self):
        '''absolute path to Darned database file'''
        return os.path.abspath(self.__darned_path)
        
    def db_name(self):
        '''Darned db name'''
        return os.path.basename(self.__darned_path)
                 
    def size(self):
        '''number of the Darned entories'''
        return self.__db_size

        
class VCFReader(object):
    '''
    VCFReader class provides that list of VCF file and utils methods
    >>> vr = VCFReader(path_to_vcf_file)
    >>> vr.db
    Returns array of vcf file
    '''
    
    def __init__(self, filename):
        self.__vcf = filename
        self.db = self.__generate_vcf_set()

    def __generate_vcf_set(self):
        _vcf_recs = []
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
        '''number of entory of the parsed vcf records'''
        return self.__size

    def vcf_name(self):
        return os.path.basename(self.__vcf)

    def editing_types(self):
        '''All type of the base substitutions'''
        return self.__substitutions
        
    def ag_count(self):
        '''A-to-G editing alone'''
        return self.__substitutions.get('A-to-G')

    def target_count(self, types):
        '''Specify type of base substitutions'''
        return self.__substitutions.get(types)
        
    def other_mutations_count(self):
        i = 0
        for k in self.__substitutions:
            if not k == 'A-to-G':
                i += self.__substitutions[k]
        return i
        

class __CSVReader(object):
    '''
    CSVReader class provides to generate array of CSV file
    >>> csv = VCFReader(path_to_csv_file)
    >>> csv.db
    Returns array of csv file
    '''
    
    def __init__(self, filename):
        self.__filename = filename
        self.db = self.__generate_csv_set()
        
    def __generate_csv_set(self):
        csv_recs = []
        with open(self.__filename) as f:
            for line in f:
                if not line.startswith("track") and not line.startswith('#'):
                    rec = line.split(',')
                    if rec[0].startswith('chr'):
                        chrom = re.sub(r'^chr', '', rec[0], 1)
                    if rec[1].find(','):
                        pos = rec[1].replace(',', '')
                    csv_recs.append(chrom + ':' + pos)
        self.__size = len(csv_recs)
        return csv_recs

    def size(self):
        return self.__size

    def name(self):
        return os.path.basename(self.__filename)

    def ag_count(self):
        raise NotImplementedError

    def other_mutations_count(self):
        raise NotImplementedError
        

#if __name__ == '__main__':
#    import sys
#    csv = __CSVReader(sys.argv[1])
#    print csv.size()
#    print csv.name()
#    print csv.db[:10]
#    


class Benchmark(object):
    '''
    >>> darned_db = DarnedReader(sp='human_hg19', source='Brain', db='Path_to_Darned_DB')
    >>> editing_db = VCFReader(filename)
    >>> bench = Benchmark(answer=darned_db, predict=candidate_db)
    >>> bench.answer
    returns set opf darned
    >>> bench.predict
    returns set of candidate sites from vcf
    >>> bench.intersect
    returns set of the intersection between sets
    '''
    
    def __init__(self, answer=None, predict=None):

        # TODO: fix here! to convert list to set
        self.answer = set([_.split(":")[:2] for _ in answer])
        self.predict = set([_.split(":")[:2] for _ in predict])
        
        if len(self.answer) == 0:
            raise ValueError, 'Answer data set has no entory'
        elif len(self.predict) == 0:
            raise ValueError, 'Candidate data set has no entory'
                
        self.intersect = self.answer.intersection(self.predict)

    def __str__(self):
        return "Answer set[%d], Candidate set[%d]\n" % (len(self.answer), len(self.predict))

    def precision(self):
        try:
            _precision = len(self.intersect)/len(self.predict)
            return _precision
        except ZeroDivisionError:
            pass
        finally:
            return 0
            
    def recall(self):
        try:
            _recall = len(self.intersect)/len(self.answer)
            return _recall
        except ZeroDivisionError:
            pass
        finally:
            return 0
            
    def f_measure(self):
        _precision = self.precision()
        _recall = self.recall()
        try:
            _f = 2*_recall*_precision/(_recall+_precision)
            return _f
        except ZeroDivisionError:
            pass
        finally:
            return 0
