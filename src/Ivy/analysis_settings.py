# Setting of `ivy`
IVY_SETTINGS =  {
    'BASIC_OPT': {
        'FASTA': None,
        'RNA_BAM': None,
        'DNA_BAM': None,
        'OUT_NAME': 'ivy.vcf',
        'REGIONS': 'All',
        'GTF': None,
        'IS_SINGLE': False,
        'ONE_BASED': False,
        'N_THREADS': 1,
        'DRY_RUN': False,
        'VERBOSE': False,
    },
    
    'SAMPLE_OPT': {
        'STRAND': False,
        'KO_STRAIN': False,
        'REPLICATE': False,
    },

    'FILTERS': {
        
        'BASIC_FILT': {
            'MIN_MUT_FREQ': 0.1,
            'MIN_RNA_COV': 10,
            'MIN_DNA_COV': 10,
            'RM_DUPLICATED': True,
            'RM_DELETION': True,
            'RM_INSERTION': True,
            'MIN_RNA_MAPQ': 30,
            'MIN_DNA_MAPQ': 28,
            'MIN_RNA_BAPQ': 28,
            'NUM_TYPE': 1,
        },
        
        'STAT_FILT': {
            'SIG_LEVEL' :0.05,
            'BAQ_BIAS' :False,
            'POS_BIAS' :False,
            'STRAND_BIAS' :False,
        },
        
        'EXT_FILT': {
            'BLAT': False,
            'SNP': False,
            'SS_NUM': False,
            'TRIM_N': False,
            'MASK_REPEAT': False,
        }
    }
}

# Setting of `edit_bench`
EDIT_BENCH_SETTINGS = {
    'APP': {
        'VCF': None,
        'CSV': None,
        'SOURCE': 'All',
        'SP': None,
        'PLOT': False,
        'OUT': 'edit_bench.log',
    },
    
    'DATA': {
        'human_hg19': 'http://darned.ucc.ie/static/downloads/hg19.txt',
        'human_hg18': 'http://darned.ucc.ie/static/downloads/hg18.txt',
        'mice_mm9': 'http://darned.ucc.ie/static/downloads/mm9.txt',
        'mice_mm10': 'http://darned.ucc.ie/static/downloads/mm10.txt',
        'fly_dm3': 'http://darned.ucc.ie/static/downloads/dm3.txt',
    }
}

