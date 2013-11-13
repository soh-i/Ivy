Ivy
===

Ivy: software package for identification of RNA editing sites based on massively parallel sequencing data

## System requirements
* Python 2.7+ (recomended)
* Unix like operating system (tested by OS X, RedHat)

## Dependencies
* Pysam version 0.7.5
* PyVCF version 0.6.4
* Matplotlib version 1.3.1

## Installation
```
python setup.py test
python setup.py install
```

## ivy
### Usage
```
ivy  -l chr21:47721030-47721057 -f reference.fa  -r rna.bam
```

#### Options:
	--version                             show program's version number and exit
	-h, --help                            show this help message and exit
	-f FASTA                              Reference genome [fasta]
	-r R_BAMS                             RNA-seq file(s) [bam]
	-b D_BAMS                             DNA-seq file(s) [bam]
	-o OUTNAME                            Output filename
	-l REGIONS                            Explore specify region [chr:start-end]
	-G GTF                                GTF file
	--one-based                           Genomic coordinate
	--num-threads=N_THREADS               Number of threads [default: 1]
	--dry-run                             Dry run ivy
	--verbose                             Show verbously messages

#### Extended filter options:
    --blat-collection                   Reduce mis-alignment with Blat [default: False]
    --snp=SNP_FILE                      Exclude variation sites [vcf]
    --ss-num=SS_NUM                     Exclude site around the splice sistes [default: 5bp]
    --trim-n=TRIM_N                     Do not call Nbp in up/down read [default: 10bp]
    --mask-repeat                       Mask repeat sequence [default: False]

####  Sample options:
    --strand                            Strand-specific seq. data is used. [default: False]
    --ko-strain                         Adar null strain is used. [default: False]
    --replicate                         Biological replicate is used [default: False]

####  Basic filter options:
    --min-ag-ratio=AG_RATIO             Min A-to-G edit base ratio [default: 0.1]
    --min-rna-coverage=MIN_RNA_COV      Min RNA read coverage [default: 10]
    --min-dna-coverage=MIN_DNA_COV      Min DNA read coverage [default: 20]
    --rm-duplicated-read=IS_DUPLICATED  Remove duplicated reads [default: True]
    --rm-deletion-read=IS_DELETION      Remove deletion reads [default: True]
    --min-mapq=MIN_MAPQ                 Min mapping quality [default: 30]
    --num-allow-type=NUM_TYPE           Number of allowing base modification type [default: 1]
    --min-baq-rna=MIN_BAQ_R             Min base call quality in RNA [default: 28]
    --min-baq-dna=MIN_BAQ_D             Min base call quality in DNA [default: 28]

####  Statistical filter options:
    --sig-level=SIG_LEVEL               Significance level [default: 0.05]
    --base-call-bias                    Consider base call bias [default: True]
    --strand-bias                       Consider strand bias [default: True]
    --positional-bias                   Consider positional bias [default: True]


## ivy_benchmark
### Benchmarking test for detection accuracy
ivy_benchmark performs a benchmarking test of detected A-to-I editing sites in _H. sapiens_, _M. musclus_, _D. melanogaster_ to evaluate accuracy using recall, recision and F-measure score, which are defined as follows: (...).

For benchmarking test, we used DARNED database (Kiran _et al_., 2013, http://beamish.ucc.ie) as true known editing sites. DARNED db  provides the comprehensive list of the previously identified RNA editing sites by researchers.

### Usage    
Called A-to-I editing sites with human brain transcriptome with hg19 genome version, following arguments:

```
ivy_benchmark --vcf test1.vcf test2.vcf --sp human_hg19 --source brain --plot
```
then, recall/precision plot outputs in your directory as pdf file.


* `--vcf`: Set VCF(variant call format) files [required]
* `--sp`: Set joined string of species and genome version [required]
* `--source`: Set specific smaple/tissues/cell line name [default: All]
* `--plot`: Use this option, to visualize precision/recall plot [default: False] 
* `--help`: Show help messages
* `--version`: Show ivy_benchmark.py version


