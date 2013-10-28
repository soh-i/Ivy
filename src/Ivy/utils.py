import os.path

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

