from setuptools import setup, find_packages

__author__ = 'Soh Ishiguro'
__email__ = 'yukke@g-language.org'
__version__ = '0.0.1'
__pkgname__ = 'Ivy'

required_packages = []
# python 2.6 does not have argsparse
try:
    import argparse
except ImportError:
    required_packages.append('argparse')

    
if __name__ == '__main__':
    setup(
        version=__version__,
        description='software package for identification of RNA editing sites based on massively parallel sequencing data',
        long_description='',
        author=__author__,
        author_email=__email__,
        url='https://github.com/soh-i/Ivy',
        name=__pkgname__,
        scripts=['bin/edit_bench', 'bin/ivy'],
        packages=find_packages('src'),
        package_dir={'': 'src'},
        entry_points = {
            'console_scripts': [
                'ivy = run_ivy:run',
                'edit_bench = run_editbench:run',
            ],
        },
        install_requires=['pysam==0.7.5', 'PyVCF==0.6.4', 'matplotlib==1.3.1'],
        test_suite='Ivy.tests',
        license='not yet',
        )
