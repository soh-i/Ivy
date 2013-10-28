from Ivy.benchmark import DarnedReader, VCFReader, Benchmark

sp = 'human_hg19'
sample = 'brain'
data = '../data/test_data.vcf'
ans = DarnedReader(sp=sp, source="brain")
vcf = VCFReader(data)
bench = Benchmark(answer=ans.db, predict=vcf.db)

print vcf.size()
print vcf.vcf_name()
print vcf.editing_types()
print vcf.ag_count()
print vcf.target_count("A-to-C")
print vcf.other_mutations_count()

