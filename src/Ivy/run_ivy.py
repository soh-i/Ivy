from Ivy.alignment.stream import RNASeqAlignmentStream, DNASeqAlignmentStream
from Ivy.commandline.parse_ivy_opts import CommandLineParser
from Ivy.annotation.writer import VCFWriteHeader
from Ivy.seq import Fasta, Timer, decode_chr_name_from_file
from Ivy.version import __version__
from pprint import pprint as p
import subprocess
import multiprocessing
import logging
import string
import shutil
import os, os.path, sys, re
import time

__program__ = 'run_ivy'
__author__ = 'Soh Ishiguro <yukke@g-language.org>'
__license__ = ''
__status__ = 'development'
__version__ = __version__

# Global configuration
logger = logging.getLogger(__name__)
parse = CommandLineParser()
params = parse.ivy_parse_options()
vcf = VCFWriteHeader(params)

def run():
    if params.n_threads == 1 or not params.n_threads:
        vcf.make_vcf_header()
        _single_run()
    elif params.n_threads > 2:
        vcf.make_vcf_header()
        _thread_run()
    logger.debug("Finished Ivy run!")
    
def _single_run():
    '''
    Working on single thread
    '''

    #vcf.make_vcf_header()
    logger.debug("Beginning Ivy run (v." + __version__ + ")" )
    
    # RNA-seq data
    if params.r_bams:
        logger.debug("Loading RNA-seq bam file '{0}'".format(params.r_bams))
        rna_pileup_alignment = RNASeqAlignmentStream(params)
        #print dir(rna_pileup_alignment)
        for rna in rna_pileup_alignment.filter_stream():
            print Printer.to_tab(rna)
            
    # DNA-Seq data
    if params.d_bams:
        logger.debug("Loading DNA-seq bam file '{0}'".format(params.d_bams))
        dna_pileup_alignment = DNASeqAlignmentStream(params)
        for dna in dna_pileup_alignment.filter_stream():
            print '{chrom:}\t{pos:}\t{ref:}\t{alt:}'.format(
                chrom=dna.get('chrom'),
                pos=dna.get('pos'),
                ref=dna.get('ref'),
                alt=dna.get('alt'))
    #with open(params.outname, 'w') as f:
    #f.write()
    logger.debug("Finished Ivy run!")

def _thread_run():
    '''
    Wrapper of working on multiprocessing
    '''
    
    logger.debug("Beginning Ivy run v." + __version__)
    save_path = './block_fasta'
    
    # create tmp directory
    if not os.path.isdir(save_path):
        os.mkdir(save_path)
    fasta_files = _get_fa_list(save_path)
    
    # generate worker
    if len(fasta_files) != params.n_threads:
        logger.debug("Remove old splited fasta")
        shutil.rmtree(save_path)
        
        # split fasta by number of threads
        logger.debug("Split reference genome by number of threads")
        fa = Fasta(fa=params.fasta)
        fa.split_by_blocks(params.n_threads)
        
        # call workers
        with Timer() as t:
            __start_worker(params.n_threads, _get_fa_list(save_path))
        
    elif len(fasta_files):
        logger.debug("Used existing splited reference genome")
        with Timer() as t:
            __start_worker(params.n_threads, _get_fa_list(save_path))
            
    logger.debug("Start to merge files")
    __merge_tmp_files(tmp_path='ivy_tmp')
    
    logger.debug("Elapsed time: {0}d {0}h {1}min {2}sec".format(
        int(t.interval*1/24*1/3600), int(t.interval*1/3600), int(t.interval*1/60), int(t.interval)))
        
def __multi_pileup(seq_files):
    '''
    Call pileup methods to run as multiprocessing
    '''
    
    current = multiprocessing.current_process()
    logger.debug("Process: {0}".format(current))
    logger.debug("PID: {0}".format(current.pid))
    logger.debug("RUN: {0}".format(current.run))
    logger.debug("Identity: {0}".format(current._identity))
    #logger.debug("Class: {0}".format(current.start.im_func.func_name))
    
    chromosome_list = decode_chr_name_from_file([seq_files])
    path_to_tmp = 'ivy_tmp'
    for seq_files in chromosome_list:
        reverted = seq_files[0] + "_" + "-".join(seq_files[1:]) + ".fa"
        params.fasta = os.path.join('block_fasta', reverted)
        logger.debug("Subjected fasta: '{0}'".format(params.fasta))
        
        for chrom in seq_files[1:]: # Skip serial number in the 1st element
            params.region.chrom = chrom
            pileup_iter = RNASeqAlignmentStream(params)
            if not os.path.isdir(path_to_tmp):
                os.mkdir(path_to_tmp)
                
            out = open(path_to_tmp+'/tmp_'+chrom+'.log', mode='w', buffering=False)
            try:
                for p in pileup_iter.filter_stream():
                    out.write(Printer.to_tab(p)+'\n')
                    
            except KeyboardInterrupt:
                logger.debug("Interrupted Ivy run, PID: {0}".format(current.pid))
                multiprocessing.terminate()
                
            finally:
                out.close()
                
    logger.debug("Finished to pileup in '{0}'".format(chrom))

def __merge_tmp_files(tmp_path=None):
    '''
    Merge temporary files into a VCF as final result
    
    Args:
     tmp_path(str): path to Ivy tmp directory
    '''
    
    if not os.path.isdir(tmp_path):
        print "Path: {0} is not found".format(tmp_path)
        return False
    
    file_list = []
    for f in os.listdir(tmp_path):
        p = re.compile(r'^tmp_{1}(.+).{1}log$')
        ma = p.match(f)
        f_name =  ma.group(1)
        if f_name.startswith('chr'):
            f_name = f_name.replace('chr', '', 1)
            if f_name.isdigit():
                file_list.append(int(f_name))
            else:
                file_list.append(str(f_name))
        else:
            if f_name.isdigit():
                file_list.append(int(f_name))
            else:
                file_list.append(str(f_name))
                
    revert_files = [tmp_path + '/tmp_chr' + str(_) + '.log' for _ in sorted(file_list)]
    tmp_size = len(revert_files)
    
    if not os.access('/bin/cat', os.X_OK):
        print "Could not merge analysis files (Not found executalble 'cat' command)"
        return False
    command = "/bin/cat "
    command += " ".join(revert_files)
    name = "ivy_merged_result.log"
    fout = open(name, mode='w', buffering=False)
    subprocess.call([command], stdout=fout, shell=True)
    logger.debug("Finished to {0} merge tmp files to {1}".format(tmp_size, name))

def __start_worker(cpus, fas):
    '''
    Start Ivy child processes for multiprocessing run
    
    Args:
     cpus(int): Number of used cpus for run
     fas(list): Fasta files get by _get_fa_list(path)
    '''
    if cpus < 1 and len(fas) < 1:
        print "Number of cpus or seq. len. is too small"
        return False
        
    logger.debug("Number of CPUs: {:d}".format(cpus))
    p = multiprocessing.Pool(cpus)
    seq = p.map(__multi_pileup, fas)
    
def list_fasta_files(path, suffix):
    '''
    Args:
     path(str): path to directory which contains fasta files
     suffix(str): suffix in fasta files, default sets '.fa'
    Return:
     fasta files(list) in given path
    '''
    
    if suffix:
        # uses custamized suffix of fasta
        return [_ for _ in os.listdir(path) if _.endswith(ext)]
    else:
        # default
        return [_ for _ in os.listdir(path) if _.endswith('.fa')]

def _get_fa_list(path):
    '''
    Args:
     path(str): blocked fasta contained directory
    Returns:
     files(list): only .fa file which located in given path
    '''
    
    fa = []
    for _ in os.listdir(path):
        if _.endswith(".fa"):
            fa.append(_)
    return fa

            
class Printer(object):
    @staticmethod
    def pprint(data, *args, **kwargs):
        '''
        Simple pretty print yielded pileuped data
        '''
        print '{'
        for key in data:
            print '  {key}({types}): {val}'.format(
                key=key, val=data[key], types=type(data[key]).__name__)
        print '}'
 
    @staticmethod
    def to_tsv(data, *args, **kwargs):
        entory = ''
        for key in data:
            if isinstance(data[key], str):
                entory += data[key]+"\t"
            elif isinstance(data[key], int) or isinstance(data[key], float):
                entory += str(data[key])+"\t"
            elif isinstance(data[key], tuple):
                entory += ','.join([str(_) for _ in data[key]]) + "\t"
        return entory
         
    @staticmethod
    def to_tab(data, *args, **kwargs):
        return '{chrom:}\t{pos:}\t{ref:}\t{alt:}\t{dp4:}\t{alf:}'.format(
            chrom=data.get('chrom'),
            pos=data.get('pos'),
            ref=data.get('ref'),
            alt=data.get('alt'),
            dp4=",".join([str(_) for _ in data.get('dp4')]),
            alf=data.get('allele_freq'))
        
         
if __name__ == '__main__':
    print "Please execute bin at APP_ROOT/bin/ivy"
