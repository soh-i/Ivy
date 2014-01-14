Ivy: software package for detection of RNA editing sites based on high-throughput sequencing data.
====

Ivy is a software package that faster and accurately detects RNA editing events in higher eukaryotes (such as human, mouse, and melanogaster) RNA-seq data.

Ivy includes 2 command line tools:

* `ivy` is detection tool for RNA editing site.
* `edit_bench` is benchmark tool to test the detection accuracy.


# System requirements
* Python 2.7+ (v2.7.5 recomended)
* Unix like operating system (tested by OS X 10.9 and RedHat)

# Dependencies
The required dependencies to install the Ivy is `pysam`, `pyVCF`, `fisher`, and `matplotlib`.

# Installation
```
pip install numpy
pip install matplotlib
python setup.py install
```

# Documentation
* [Ivy v0.0.1: documentation and manual](http://web.sfc.keio.ac.jp/~t10078si/ivy.html)

# Licence
Ivy is available under the GPLv2 license. See the `LICENSE` file for more information.