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
### `ivy_benchmark.py`: test for detection accuracy
`ivy_benchmark.py` performs a benchmarking test of detected A-to-I editing sites in _H. sapiens_, _M. musclus_, _D. melanogaster_ to evaluate accuracy using recall, recision and F-measure score. We used DARNED database (Kiran _et al_., 2013) as true known editing sites, DARNED provides the comprehensive list of the previously described RNA editing sites.
    
Called A-to-I editing sites in human brain sample with hg19 genome version, following arguments:

```
ivy_benchmark.py --vcf test.vcf --sp human_hg19 --source brain
```
* `--vcf`: Set VCF(variant call format) file [required]
* `--sp`: Set joined string of species and genome version [required]
* `--source`: Set specific smaple/tissues/cell line name [default: All]
* `--plot`: Use this option, to visualize precision/recall plot [default: False] 
* `--help`: Show help messages
* `--version`: Show ivy_benchmark.py version

## Dependencies
* `pip install -r requirements.txt` resoluves all dependencies.
* Pysam version 0.7.5
* PyVCF version 0.6.4
* Numpy version 1.7.1
* Matplotlib version 1.3.1
