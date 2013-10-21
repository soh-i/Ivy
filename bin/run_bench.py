from Ivy.benchmark import *
import argparse

if __name__ == '__main__':
    version = '0.0.1'
    desc = "Benchmarking test for detected RNA editing sites"
    parser = argparse.ArgumentParser(version=version, description=desc)
    parser.add_argument('--vcf',
                        required=True,
                        dest='vcf_file',
                        action='store',
                        type=str,
                        help='set VCF file'
                    )
    args = parser.parse_args()
    
    if args.vcf_file:
        v = VCFReader(args.vcf_file)
        print v.vcf_name()
        print v.cnt()
        print v.editing_types()
        print v.ag_count()
        print v.other_mutations_count()
