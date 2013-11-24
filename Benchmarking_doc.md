EditBench
===

EditBench: Simple python program for benchmarking test against identified RNA editing sites based on high-throughput sequencing data.

## Requirements
* Python 2.7+ (recomended)
* Unix like operating system

## Dependencies
* PyVCF version 0.6.4
* Matplotlib version 1.3.1

## Installation
```
python setup.py install
```

## Usage
ivy_benchmark performs a benchmarking test of detected A-to-I editing sites in _H. sapiens_, _M. musclus_, _D. melanogaster_ to evaluate accuracy using recall, recision and F-measure score, which are defined as follows:   

	Recall = 
	Precision = 
	F-measure = 

For benchmarking test, we used DARNED database (Kiran _et al_., 2013, http://beamish.ucc.ie) as true known editing sites. DARNED db  provides the comprehensive list of the previously identified RNA editing sites by researchers.

Called A-to-I editing sites with human brain transcriptome with hg19 genome version, following arguments:

```
ivy_benchmark --vcf test1.vcf test2.vcf --sp human_hg19 --source brain --plot
```
then, recall/precision plot outputs in your directory as pdf file.

### Required argumentss
* `--vcf`: VCF(variant call format) files [required]
* `--sp`: Species and genome version [required]

### Options
* `--source`: Specific smaple/tissues/cell line name [default: All]
* `--plot`: Plot precision/recall plot [default: False]
* `--out`: Outpust filename. [default: stdout]
* `--help`: Show help messages
* `--version`: Show ivy_benchmark.py version

