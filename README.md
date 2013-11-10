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

## Benchmarking test
### ivy_benchmark: detection accuracy testing
ivy_benchmark performs a benchmarking test of detected A-to-I editing sites in _H. sapiens_, _M. musclus_, _D. melanogaster_ to evaluate accuracy using recall, recision and F-measure score, which are defined as follows: (...).

For benchmarking test, we used DARNED database (Kiran _et al_., 2013, http://beamish.ucc.ie) as true known editing sites. DARNED db  provides the comprehensive list of the previously identified RNA editing sites by researchers.

#### Usage    
Called A-to-I editing sites with human brain transcriptome with hg19 genome version, following arguments:

```
ivy_benchmark --vcf test.vcf --sp human_hg19 --source brain --plot
```

Recall/Precision plot outputs in your directory as PDF.
![](https://f.cloud.github.com/assets/1855860/1508063/08482826-49ce-11e3-81a8-ca52e1eec73a.png)


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
