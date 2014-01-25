from setuptools import  find_packages, setup
import sys
import os.path
from Ivy.version import __version__

__author__ = 'Soh Ishiguro'
__email__ = 'yukke@g-language.org'
__pkgname__ = 'Ivy'


# DO NOT WORK PYTHON 3
if not sys.version_info[0] == 2:
    print "Sorry, Python 3 is not support yet"
    sys.exit()

def deps():
    if os.path.isfile('requirements.txt'):
        deps =  open('requirements.txt').read().splitlines()
        return deps
    else:
        sys.exit('Not found \'requirements.txt\' file!')

        
if __name__ == '__main__':
    setup(
        version=__version__,
        description='software package for identification of RNA editing sites based on massively parallel sequencing data',
        author=__author__,
        author_email=__email__,
        url='https://github.com/soh-i/Ivy',
        name=__pkgname__,
        packages=find_packages('src'),
        package_dir={'': 'src'},
        entry_points = {
            'console_scripts': [
                'ivy = Ivy.run_ivy:run',
                'edit_bench = Ivy.run_editbench:run',
            ],
        },
        install_requires=deps(),
        test_suite='Ivy.tests',
        dependency_links = ["samtools.sourceforge.net"],
        license='GNU General Public License v2.0',
        keywords = ['Bioinformatics',
                    'RNA editing',
                    'RNA-Seq',
                    'NGS',
                    'High-throughput sequencing',
                ],
        classifiers = [
            'Development Status :: 3 - Alpha',
            'Intended Audience :: Science/Research',
            'Programming Language :: Python',
            'Topic :: Scientific/Engineering :: Bio-Informatics'
            'License :: OSI Approved :: GNU General Public License v2 (GPLv2)'
        ],
        platforms=['Linux','Unix','MacOS'],
        long_description=open('README.md').read(),
    )
    
