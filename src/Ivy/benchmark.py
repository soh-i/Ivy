from __future__ import division
import vcf
import os.path
import re
import ConfigParser
from collections import Counter
import utils

class DarnedDataGenerator(object):
    '''
    DarnedDataGenerator provides to prepare data that are used for benchmarking test.
    '''
    
    def __init__(self, species=None):
        self.species = species

    def fetch_darned(self):
        '''
        Fetch specify raw dataest from darned.ucc.ie/static/downloads/ into APP_ROOT/data,
        species name must be given, and acceptable type is defined as:
        human_hg18/hg19, mice_mm9/mm10, fly_dm3
        '''
        
        __species = {
            'human_hg19':'http://darned.ucc.ie/static/downloads/hg19.txt',
            'human_hg18':'http://darned.ucc.ie/static/downloads/hg18.txt',
            'mice_mm9':'http://darned.ucc.ie/static/downloads/mm9.txt',
            'mice_mm10':'http://darned.ucc.ie/static/downloads/mm10.txt',
            'fly_dm3':'http://darned.ucc.ie/static/downloads/dm3.txt'
        }
        
        self.filename = self.species + '.txt'
        if os.path.isfile(self.filename):
            return False
            
        try:
            url = __species[self.species]
        except:
            raise RuntimeError, "Given [%s] is not valid species name" % (self.species)
        
        req = Request(url)
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
            root_path = find_app_root()
            if not os.path.isdir(root_path + '/data'):
                os.makedirs(root_path + '/data')
                print "Make directories [%s]" % (root_path + '/data')
                
            print "Dowloading [%s] from [%s] ..." % (self.filename, url)
            with open(root_path+ '/data'+ self.filename, "w") as fout:
                fout.write(response.read())
            return True

    def darned_to_csv(self, filename):
        '''
        Converting darned raw datafile to csv,
        the data that fetched from darned.ucc.ie/static/downloads/*.txt is given.
        >>> path_to_data = hg19.txt
        >>> darned_to_csv(path_to_data)
        Generate csv file into the APP_ROOT/data
        '''
        
        if not os.path.isfile(filename):
            raise RuntimeError, 'filename->[%s] is not found' % (filename)
        
        data_path = find_app_root() + '/data/'
        if not os.path.isdir(data_path):
            print "Create data dir"
            os.makedirs(data_path)
        
        name, ext = os.path.splitext(filename)
        out_name = data_path + name + '.csv'
        if os.path.isfile(out_name):
            print "File is already exisit"
            return False
        
        reader = csv.reader(open(filename, 'r'), delimiter="\t", quotechar="|")
        try:
            line_n = 0
            out = open(out_name, 'w')
            for row in reader:
                line_n += 1
                source = row[8]
                if len(source):
                    mod = source.replace(r';', ',').replace(r',', ';').replace(r'; ',';').replace(r' ','_').replace(r'_T','T')
                    out.write(",".join(row[:8]) + ",")
                    out.write(mod + ",")
                    out.write(",".join(row[9:]) + "\n")
        except:
            raise ValueError, ('Parsing error at line No.[%d]') % (line_n)
        
        finally:
            out.close()

class DarnedReader(object):
    '''
    DarnedReader generates the subset of DARNED db
    >>> db = DarnedReader(sp='human_hg19', source='Brain', db='Path_to_Darned_DB')
    Returns list of subset of darned db
    '''
    
    def __init__(self, sp='', source=None, db_path=None):
        self.__sp = sp
        if source is not None:
            self.source = source
        else:
            self.source = None
        """
        if db_path is None:
            conf_path = utils.find_app_root()
            conf = ConfigParser.RawConfigParser()
            conf.read(conf_path+ '/data.ini')

            if conf.has_section('DARNED'):
                sp = conf.get('DARNED', self.__sp)
                self.__darned = {self.__sp: conf_path+ sp}
                self.db = self.__generate_darned_set()
            else:
                raise (ConfigParser.NoSectionError,
                       'Invalid species name [%s] is given' % (self.__sp))
        
        elif db_path is not None:
            self.__darnd = {self.__sp: db_path}
            self.db = self.__generate_darned_set()
        """
        if db_path is not None:
            self.__darned = {self.__sp: db_path}
        else:
            raise RuntimeError, ('db_path is not set')

    def __str__(self):
        return 'Darned file [%s]' % (self.__darned[self.__sp])
        
    def __generate_darned_set(self):
        
        # Store selected records
        self.__size = 0
        if self.__sp and self.source is not None:
            darned_list = []
            with open(self.__darned[self.__sp], 'r') as f:
                for line in f:
                    if not line.startswith('chrom'):
                        data = line.split(",")
                        chrom = data[0]
                        pos = data[1]
                        darned_source = data[8]
                        
                        if darned_source == self.source:
                            darned_list.append(chrom+ ':'+ pos+ self.source)
                            self.__size += 1
                return darned_list
                
        # Store all Darned records (default)
        elif self.__sp and self.source is None \
             or self.source is 'all' or self.source is 'All':
            selected = []
            with open(self.__darned[self.__sp], 'r') as f:
                for line in f:
                    if not line.startswith('chrom'):
                        data = line.split(',')
                        chrom = data[0]
                        pos = data[1]
                        selected.append(chrom+ ':'+ pos+ 'All')
                        self.__size += 1
                return selected
                
        elif not self.__sp:
            raise (RuntimeError, 'Given species name[%s] is not valid' % (self.__sp))

    def sp(self):
        ''' given species name '''
        return self.__sp

    def path(self):
        ''' absolute path to Darned database file'''
        return os.path.abspath(self.__darned[self.__sp])
        
    def db_name(self):
        '''Darned db name'''
        return os.path.basename(self.path())
                 
    def size(self):
        '''number of the Darned entories'''
        return self.__size

        
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

    def ag_enrichment(self):
        raise NotImplementedError

        
