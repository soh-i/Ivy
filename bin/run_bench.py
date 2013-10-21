from Ivy.benchmark import *
import argparse

if __name__ == '__main__':
    version = '0.0.1'
    desc = "Benchmarking test for detected RNA editing sites"
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument('--vcf',
                        required=True,
                        dest='vcf_file',
                        action='store',
                        type=file,
                        help='set VCF file'
                    )
    parser.add_argument('--db',
                        required=False,
                        dest='db_file',
                        action='store',
                        type=str,
                        help='set db name'
                    )
    parser.add_argument('--config',
                        required=False,
                        dest='conf_file',
                        action='store',
                        type=file,
                        help='set config file'
                    )
    parser.add_argument('--verbose',
                        required=False,
                        dest='is_verbose',
                        action='store_true',
                        help='show verbosely output'
                    )
    parser.add_argument('--version', action='version', version=version)
    
    args = parser.parse_args()
    
    if args.vcf_file:
        v = VCFReader(args.vcf_file)
        print v.vcf_name()
        print v.cnt()
        print v.editing_types()
        print v.ag_count()
        print v.other_mutations_count()
