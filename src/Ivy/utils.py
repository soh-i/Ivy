import os.path
import csv
import re
from urllib2 import Request, urlopen, URLError

def find_app_root():
    '''
    Absolute path to your project root from setup.py location
    '''
    root = os.path.dirname(__file__)
    while not os.path.exists(os.path.join(root, 'setup.py')):
        root = os.path.abspath(os.path.join(root, os.path.pardir))
    return root

def __end_url_basename(p):
    """Returns the final component of a pathname"""
    i = p.rfind('/') + 1
    return p[i:]

def setup_darned():
    for sp in ['human_hg19', 'human_hg18', 'mice_mm9', 'mice_mm10', 'fly']:
        fetch_darned(sp)
    
def fetch_darned(species):
    '''
    Fetch specify raw dataest from darned.ucc.ie/static/downloads/ into APP_ROOT/data
    species name must be given,
    species type is defined as: human_hg18/hg19, mice_mm9/mm10, fly_dm3
    '''
    __species = {
        'human_hg19':'http://darned.ucc.ie/static/downloads/hg19.txt',
        'human_hg18':'http://darned.ucc.ie/static/downloads/hg18.txt',
        'mice_mm9':'http://darned.ucc.ie/static/downloads/mm9.txt',
        'mice_mm10':'http://darned.ucc.ie/static/downloads/mm10.txt',
        'fly_dm3':'http://darned.ucc.ie/static/downloads/dm3.txt'
    }
    
    filename = species + '.txt'
    if os.path.isfile(filename):
        return False
        
    try:
        url = __species[species]
    except:
        raise RuntimeError, "Given [%s] is not valid species name" % (species)
        
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
                
        print "Dowloading [%s] from [%s] ..." % (filename, url)
        with open(root_path+ '/data'+ filename, "w") as fout:
            fout.write(response.read())
    return True

def darned_to_csv(filename):
    '''
    Converting darned row datafile to csv,
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
        
if __name__ == '__main__':
    fetch_darned("human_hg188")
