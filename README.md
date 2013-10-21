Ivy
===

Ivy: Tools for identifying RNA editing sites based on HTSeq data

## System requirements
* Python 2.7+ (recomended)
* Unix like operating system (OS X, Redhad, and more)

## Installation
```
python setup.py test
python setup.py install
```

## Usage
### Benchmarking test
```
python bin/run_bench.py --vcf test --sp human_hg19
```

Shows other options

```
python bin/run_bench.py --help
```

## Dependencies
* Pysam version 0.7.5
* PyVCF version 0.6.4
