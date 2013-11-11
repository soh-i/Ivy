from setuptools import setup, find_packages

__author__ = 'Soh Ishiguro'
__email__ = 'yukke@g-language.org'
__version__ = '0.0.1'
__pkgname__ = 'Ivy'

setup(
    version=__version__,
    description='software package for identification of RNA editing sites based on massively parallel sequencing data',
    author=__author__,
    author_email=__email__
    url='https://github.com/soh-i/Ivy',
    name=__pkgname__,
    license='not yet',
    package=find_packages('src'),
    package_dir={'': 'src'},
    install_requires=['pysam==0.7.5', 'PyVCF==0.6.4', 'matplotlib==1.3.1'],
    test_suite='Ivy.tests',
    )
