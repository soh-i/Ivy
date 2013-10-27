import os.path
import csv
import re
from urllib2 import Request, urlopen, URLError

def find_app_root():
    root = os.path.dirname(__file__)
    while not os.path.exists(os.path.join(root, 'setup.py')):
        root = os.path.abspath(os.path.join(root, os.path.pardir))
    return root

def fetch_darned():
    if os.path.isfile(filename):
        return False

    url = 'http://beamish.ucc.ie/data_hg19.txt'
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
            root_path = os.path.dirname(filename)
            if not os.path.isdir(root_path):
                os.makedirs(root_path)
                print "Made directories [%s]" % (root_path)
                print "Dowloading [%s] from [%s] ..." % (filename, url)
                with open(filename, "w") as fout:
                    fout.write(response.read())
                return True

def darned_to_csv(filename):
    if not os.path.isfile(filename):
        raise RuntimeError, 'filename->[%s] is not found' % (filename)
        
    data_path = find_app_root() + '/data/'
    if not os.path.isdir(data_path):
        os.makedirs(data_path)
        
    name, ext = os.path.splitext(filename)
    out_name = data_path + name + '.csv'
    if os.path.isfile(out_name):
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
                out.write(",".join(row[:8]))
                out.write(mod + ",")
                out.write(",".join(row[9:]) + "\n")
    except:
        raise ValueError, ('Parsing error at line No.[%d]') % (line_n)
        
    finally:
        out.close()
        
if __name__ == '__main__':
    print darned_to_csv("hg19.txt")
