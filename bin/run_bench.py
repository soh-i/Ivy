from Ivy.benchmark import *
from optparse import OptionParser

if __name__ == '__main__':
    version = '0.0.1'
    usage = '%prog -v test.vcf'
    parser = OptionParser(usage=usage, version=version)
    parser.add_option(
        '-v',
        action='store', type='string', dest='vcf_file', help='Set vcf filename'
    )
    options, args = parser.parse_args()
    if options.vcf_file:
        v = VCFReader(options.vcf_file)
        print v.vcf_name()
        print v.cnt()
        print v.editing_types()
        print v.ag_count()
        print v.other_mutations_count()
    else:
        parser.error("VCF file must be given")

