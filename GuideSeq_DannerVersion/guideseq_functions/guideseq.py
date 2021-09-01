#!/usr/bin/env python2.7
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 15 12:53:02 2020

@author: edanner

Starting point was the Uditas software. Many functions have been tweaked and many new ones added.
The indel analysis portions are untouched



"""
from __future__ import print_function
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pylab
import pandas as pd
import numpy as np
import os
import sys
import gzip
import itertools
import operator
import subprocess
import twobitreader
from Bio.Alphabet import IUPAC
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import pysam

#not sure if I need these
class Error(Exception):
    """Base class for exceptions in this module."""
    pass

class StrandError(Error):
    """Exception raised for errors in the strand information.
    Attributes:
        expression -- input expression in which the error occurred
        message -- explanation of the error
    """
    def __init__(self, message):
        self.message = message

class ReactionTypeError(Error):
    """Exception raised for errors in the reaction type to be processed.
    Attributes:
        expression -- input expression in which the error occurred
        message -- explanation of the error
    """
    def __init__(self, message):
        self.message = message

#these are copied and unchanged from the Uditas v1 software

################################################################################
# Open .fastq or .fastq.gz files for reading
################################################################################
def open_fastq_or_gz(filename):
    if filename.endswith(".fastq") and os.access(filename, os.F_OK):
        return open(filename, "rU")
    elif filename.endswith(".fastq.gz") and os.access(filename, os.F_OK):
        return gzip.open(filename, "rb")
    elif filename.endswith(".fastq") and os.access(filename + ".gz", os.F_OK):
        return gzip.open(filename + ".gz", "rb")
    elif filename.endswith(".fastq.gz") and os.access(filename[:-3], os.F_OK):
        return open(filename[:-3], "rU")
    raise IOError("Unknown file: " + filename)

################################################################################
# Hamming distance
# From http://code.activestate.com/recipes/499304-hamming-distance/
################################################################################
def hamm_dist(str1, str2):
    assert len(str1) == len(str2)
    ne = operator.ne
    return sum(itertools.imap(ne, str1, str2))



################################################################################
# Mask sequence by quality score
################################################################################
def mask(seq, qual, min_qual=12):

    return "".join((b if (ord(q) - 33) >= min_qual else "N") for b, q in itertools.izip(seq, qual))


################################################################################
# get the reverse-complement DNA sequence
################################################################################
def reverse_complement(seq):
    seq_dict = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'N': 'N', 'a': 't', 't': 'a', 'g': 'c', 'c': 'g'}
    return "".join([seq_dict[base] for base in reversed(seq)])





################################################################################
# create list of output files
#I added a bit to allow for the function that makes a fastq of the correctly primed targets
#I  need to add to make  mishas umi file
################################################################################
    
def create_filename(dir_sample, N7, N5, filetype):
    main_folder = os.path.join(dir_sample, N7 + '_' + N5)
    if filetype == 'mainfolder':
        return main_folder
    elif filetype == 'amplicons':
        return os.path.join(main_folder, 'amplicons')
    elif filetype == 'bam_global':
        return os.path.join(main_folder, 'bam_genome_global_files')
    elif filetype == 'bam_local':
        return os.path.join(main_folder, 'bam_genome_local_files')
    
    elif filetype == 'R1fastq':
        return os.path.join(main_folder, 'fastq_files', N7 + '_' + N5 + '_R1.fastq')
    elif filetype == 'R1fastqgz':
        return os.path.join(main_folder, 'fastq_files', N7 + '_' + N5 + '_R1.fastq.gz')
    elif filetype == 'R2fastq':
        return os.path.join(main_folder, 'fastq_files', N7 + '_' + N5 + '_R2.fastq')
    elif filetype == 'R2fastqgz':
        return os.path.join(main_folder, 'fastq_files', N7 + '_' + N5 + '_R2.fastq.gz')
      
    #### I added these as the sequences that were correctly primed
    elif filetype == 'R1fastq_CorrPrime':
        return os.path.join(main_folder, 'fastq_files', N7 + '_' + N5 + '_R1.CorrPrime.fastq')
    elif filetype == 'R1fastqgz_CorrPrime':
        return os.path.join(main_folder, 'fastq_files', N7 + '_' + N5 + '_R1.CorrPrime.fastq.gz')
    elif filetype == 'R2fastq_CorrPrime':
        return os.path.join(main_folder, 'fastq_files', N7 + '_' + N5 + '_R2.CorrPrime.fastq')
    elif filetype == 'R2fastqgz_CorrPrime':
        return os.path.join(main_folder, 'fastq_files', N7 + '_' + N5 + '_R2.CorrPrime.fastq.gz')
    elif filetype == 'localpriming_data':
        return os.path.join(main_folder, 'fastq_files', N7 + '_' + N5 + 'priming.xlsx')
    elif filetype == 'priming_summary':
        return os.path.join(main_folder, 'results', 'priming_summary.xlsx')
    #####
   
    
    elif filetype == 'R1trimmed':
        return os.path.join(main_folder, 'cutadapt_files', N7 + '_' + N5 + '_R1.trimmed.fastq.gz')
    elif filetype == 'R2trimmed':
        return os.path.join(main_folder, 'cutadapt_files', N7 + '_' + N5 + '_R2.trimmed.fastq.gz')
    elif filetype == 'trimmed_report':
        return os.path.join(main_folder, 'cutadapt_files', N7 + '_' + N5 + '.trimmed.report.txt')
    elif filetype == 'sam_genome_local':
        return os.path.join(main_folder, 'sam_genome_local_files', N7 + '_' + N5 + '.sam')
    elif filetype == 'sam_report_genome_local':
        return os.path.join(main_folder, 'sam_genome_local_files', N7 + '_' + N5 + '.sam.report.txt')
    elif filetype == 'bam_genome_local':
        return os.path.join(main_folder, 'bam_genome_local_files', N7 + '_' + N5 + '.bam')
    elif filetype == 'sorted_bam_genome_local':
        return os.path.join(main_folder, 'bam_genome_local_files', N7 + '_' + N5 + '.sorted.bam')
    elif filetype == 'sorted_bai_genome_local':
        return os.path.join(main_folder, 'bam_genome_local_files', N7 + '_' + N5 + '.sorted.bam.bai')
    elif filetype == 'sam_plasmid_local':
        return os.path.join(main_folder, 'sam_plasmid_local_files', N7 + '_' + N5 + '.sam')
    elif filetype == 'sam_report_plasmid_local':
        return os.path.join(main_folder, 'sam_plasmid_local_files', N7 + '_' + N5 + '.sam.report.txt')
    elif filetype == 'bam_plasmid_local':
        return os.path.join(main_folder, 'bam_plasmid_local_files', N7 + '_' + N5 + '.bam')
    elif filetype == 'sorted_bam_plasmid_local':
        return os.path.join(main_folder, 'bam_plasmid_local_files', N7 + '_' + N5 + '.sorted.bam')
    elif filetype == 'sorted_bai_plasmid_local':
        return os.path.join(main_folder, 'bam_plasmid_local_files', N7 + '_' + N5 + '.sorted.bam.bai')
    elif filetype == 'unmapped_bam_plasmid_local':
        return os.path.join(main_folder, 'bam_plasmid_local_files', N7 + '_' + N5 + '_unmapped.bam')
    elif filetype == 'qsorted_unmapped_bam_plasmid_local':
        return os.path.join(main_folder, 'bam_plasmid_local_files', N7 + '_' + N5 + '_qsorted_unmapped.bam')
    elif filetype == 'unmapped_plasmid_R1fastq':
        return os.path.join(main_folder, 'plasmid_unmapped_fastq_files', N7 + '_' + N5 + '_plasmid_unmapped_R1.fastq')
    elif filetype == 'unmapped_plasmid_R2fastq':
        return os.path.join(main_folder, 'plasmid_unmapped_fastq_files', N7 + '_' + N5 + '_plasmid_unmapped_R2.fastq')
    elif filetype == 'unmapped_plasmid_R1fastqgz':
        return os.path.join(main_folder, 'plasmid_unmapped_fastq_files', N7 + '_' + N5 + '_plasmid_unmapped_R1.fastq.gz')
    elif filetype == 'unmapped_plasmid_R2fastqgz':
        return os.path.join(main_folder, 'plasmid_unmapped_fastq_files', N7 + '_' + N5 + '_plasmid_unmapped_R2.fastq.gz')
    elif filetype == 'sam_amplicons':
        return os.path.join(main_folder, 'sam_amplicon_files', N7 + '_' + N5 + '.sam')
    elif filetype == 'sam_report_amplicons':
        return os.path.join(main_folder, 'sam_amplicon_files', N7 + '_' + N5 + '.sam.report.txt')
    elif filetype == 'bam_amplicons':
        return os.path.join(main_folder, 'bam_amplicon_files', N7 + '_' + N5 + '.bam')
    elif filetype == 'sorted_bam_amplicons':
        return os.path.join(main_folder, 'bam_amplicon_files', N7 + '_' + N5 + '.sorted.bam')
    elif filetype == 'sorted_bai_amplicons':
        return os.path.join(main_folder, 'bam_amplicon_files', N7 + '_' + N5 + '.sorted.bam.bai')
    elif filetype == 'unmapped_bam_amplicons':
        return os.path.join(main_folder, 'bam_amplicon_files', N7 + '_' + N5 + '_amplicons_unmapped.bam')
    elif filetype == 'qsorted_unmapped_bam_amplicons':
        return os.path.join(main_folder, 'bam_amplicon_files', N7 + '_' + N5 + '_qsorted_amplicons_unmapped.bam')
    elif filetype == 'unmapped_amplicons_R1fastq':
        return os.path.join(main_folder, 'amplicons_unmapped_fastq_files', N7 + '_' + N5 + '_amplicons_unmapped_R1.fastq')
    elif filetype == 'unmapped_amplicons_R2fastq':
        return os.path.join(main_folder, 'amplicons_unmapped_fastq_files', N7 + '_' + N5 + '_amplicons_unmapped_R2.fastq')
    elif filetype == 'unmapped_amplicons_R1fastqgz':
        return os.path.join(main_folder, 'amplicons_unmapped_fastq_files',
                            N7 + '_' + N5 + '_amplicons_unmapped_R1.fastq.gz')
    elif filetype == 'unmapped_amplicons_R2fastqgz':
        return os.path.join(main_folder, 'amplicons_unmapped_fastq_files',
                            N7 + '_' + N5 + '_amplicons_unmapped_R2.fastq.gz')
    elif filetype == 'unmapped_amplicons_report':
        return os.path.join(main_folder, 'amplicons_unmapped_fastq_files', N7 + '_' + N5 + '.unmapped.report.txt')
    elif filetype == 'sam_genome_global':
        return os.path.join(main_folder, 'sam_genome_global_files', N7 + '_' + N5 + '.sam')
    elif filetype == 'sam_report_genome_global':
        return os.path.join(main_folder, 'sam_genome_global_files', N7 + '_' + N5 + '.sam.report.txt')
    elif filetype == 'bam_genome_global':
        return os.path.join(main_folder, 'bam_genome_global_files', N7 + '_' + N5 + '.bam')
    elif filetype == 'sorted_bam_genome_global':
        return os.path.join(main_folder, 'bam_genome_global_files', N7 + '_' + N5 + '.sorted.bam')
    elif filetype == 'sorted_bai_genome_global':
        return os.path.join(main_folder, 'bam_genome_global_files', N7 + '_' + N5 + '.sorted.bam.bai')
    elif filetype == 'results_amplicons':
        return os.path.join(main_folder, 'results', N7 + '_' + N5)  # We will append the window size later
    elif filetype == 'results_plasmid':
        return os.path.join(main_folder, 'results', N7 + '_' + N5 + '_results_plasmid.xlsx')
    elif filetype == 'results_all_amplicons':
        return os.path.join(main_folder, 'results', N7 + '_' + N5 + '_results_all_amplicons.xlsx')
    elif filetype == 'results_genomewide':
        return os.path.join(main_folder, 'results', N7 + '_' + N5 + '_results_genomewide.xlsx')
    elif filetype == 'summary_all_alignments':
        return os.path.join(main_folder, 'results', N7 + '_' + N5 + '_summary_all_alignments.xlsx')
    elif filetype == 'read_counts':
        return os.path.join(main_folder, 'results', N7 + '_' + N5 + '_read_counts.xlsx')
    ## more added
    elif filetype == 'results_pipeline2_global':
        return os.path.join(main_folder, 'results', N7 + '_' + N5 + '_results_pipeline2_global.xlsx')
    elif filetype == 'results_pipeline2_local':
        return os.path.join(main_folder, 'results', N7 + '_' + N5 + '_results_pipeline2_local.xlsx')
    elif filetype == 'R1trimmed':
        return os.path.join(main_folder, 'cutadapt_files', N7 + '_' + N5 + '_R1.trimmed.fastq.gz')
    elif filetype == 'R2trimmed':
        return os.path.join(main_folder, 'cutadapt_files', N7 + '_' + N5 + '_R2.trimmed.fastq.gz')
    elif filetype == 'trimmed_report':
        return os.path.join(main_folder, 'cutadapt_files', N7 + '_' + N5 + '.trimmed.report.txt')
    elif filetype == 'sam_genome_local':
        return os.path.join(main_folder, 'sam_genome_local_files', N7 + '_' + N5 + '.sam')
    elif filetype == 'sam_report_genome_local':
        return os.path.join(main_folder, 'sam_genome_local_files', N7 + '_' + N5 + '.sam.report.txt')
    elif filetype == 'bam_genome_local':
        return os.path.join(main_folder, 'bam_genome_local_files', N7 + '_' + N5 + '.bam')
    elif filetype == 'sorted_bam_genome_local':
        return os.path.join(main_folder, 'bam_genome_local_files', N7 + '_' + N5 + '.sorted.bam')
    elif filetype == 'sorted_bai_genome_local':
        return os.path.join(main_folder, 'bam_genome_local_files', N7 + '_' + N5 + '.sorted.bam.bai')
    elif filetype == 'sam_plasmid_local':
        return os.path.join(main_folder, 'sam_plasmid_local_files', N7 + '_' + N5 + '.sam')
    elif filetype == 'sam_report_plasmid_local':
        return os.path.join(main_folder, 'sam_plasmid_local_files', N7 + '_' + N5 + '.sam.report.txt')
    elif filetype == 'bam_plasmid_local':
        return os.path.join(main_folder, 'bam_plasmid_local_files', N7 + '_' + N5 + '.bam')
    elif filetype == 'sorted_bam_plasmid_local':
        return os.path.join(main_folder, 'bam_plasmid_local_files', N7 + '_' + N5 + '.sorted.bam')
    elif filetype == 'sorted_bai_plasmid_local':
        return os.path.join(main_folder, 'bam_plasmid_local_files', N7 + '_' + N5 + '.sorted.bam.bai')
    elif filetype == 'unmapped_bam_plasmid_local':
        return os.path.join(main_folder, 'bam_plasmid_local_files', N7 + '_' + N5 + '_unmapped.bam')
    elif filetype == 'qsorted_unmapped_bam_plasmid_local':
        return os.path.join(main_folder, 'bam_plasmid_local_files', N7 + '_' + N5 + '_qsorted_unmapped.bam')
    elif filetype == 'unmapped_plasmid_R1fastq':
        return os.path.join(main_folder, 'plasmid_unmapped_fastq_files', N7 + '_' + N5 + '_plasmid_unmapped_R1.fastq')
    elif filetype == 'unmapped_plasmid_R2fastq':
        return os.path.join(main_folder, 'plasmid_unmapped_fastq_files', N7 + '_' + N5 + '_plasmid_unmapped_R2.fastq')
    elif filetype == 'unmapped_plasmid_R1fastqgz':
        return os.path.join(main_folder, 'plasmid_unmapped_fastq_files', N7 + '_' + N5 + '_plasmid_unmapped_R1.fastq.gz')
    elif filetype == 'unmapped_plasmid_R2fastqgz':
        return os.path.join(main_folder, 'plasmid_unmapped_fastq_files', N7 + '_' + N5 + '_plasmid_unmapped_R2.fastq.gz')
    elif filetype == 'sam_amplicons':
        return os.path.join(main_folder, 'sam_amplicon_files', N7 + '_' + N5 + '.sam')
    elif filetype == 'sam_report_amplicons':
        return os.path.join(main_folder, 'sam_amplicon_files', N7 + '_' + N5 + '.sam.report.txt')
    elif filetype == 'bam_amplicons':
        return os.path.join(main_folder, 'bam_amplicon_files', N7 + '_' + N5 + '.bam')
    elif filetype == 'sorted_bam_amplicons':
        return os.path.join(main_folder, 'bam_amplicon_files', N7 + '_' + N5 + '.sorted.bam')
    elif filetype == 'sorted_bai_amplicons':
        return os.path.join(main_folder, 'bam_amplicon_files', N7 + '_' + N5 + '.sorted.bam.bai')
    elif filetype == 'unmapped_bam_amplicons':
        return os.path.join(main_folder, 'bam_amplicon_files', N7 + '_' + N5 + '_amplicons_unmapped.bam')
    elif filetype == 'qsorted_unmapped_bam_amplicons':
        return os.path.join(main_folder, 'bam_amplicon_files', N7 + '_' + N5 + '_qsorted_amplicons_unmapped.bam')
    elif filetype == 'unmapped_amplicons_R1fastq':
        return os.path.join(main_folder, 'amplicons_unmapped_fastq_files', N7 + '_' + N5 + '_amplicons_unmapped_R1.fastq')
    elif filetype == 'unmapped_amplicons_R2fastq':
        return os.path.join(main_folder, 'amplicons_unmapped_fastq_files', N7 + '_' + N5 + '_amplicons_unmapped_R2.fastq')
    elif filetype == 'unmapped_amplicons_R1fastqgz':
        return os.path.join(main_folder, 'amplicons_unmapped_fastq_files',
                            N7 + '_' + N5 + '_amplicons_unmapped_R1.fastq.gz')
    elif filetype == 'unmapped_amplicons_R2fastqgz':
        return os.path.join(main_folder, 'amplicons_unmapped_fastq_files',
                            N7 + '_' + N5 + '_amplicons_unmapped_R2.fastq.gz')
    elif filetype == 'unmapped_amplicons_report':
        return os.path.join(main_folder, 'amplicons_unmapped_fastq_files', N7 + '_' + N5 + '.unmapped.report.txt')
    elif filetype == 'sam_genome_global':
        return os.path.join(main_folder, 'sam_genome_global_files', N7 + '_' + N5 + '.sam')
    elif filetype == 'sam_report_genome_global':
        return os.path.join(main_folder, 'sam_genome_global_files', N7 + '_' + N5 + '.sam.report.txt')
    elif filetype == 'bam_genome_global':
        return os.path.join(main_folder, 'bam_genome_global_files', N7 + '_' + N5 + '.bam')
    elif filetype == 'sorted_bam_genome_global':
        return os.path.join(main_folder, 'bam_genome_global_files', N7 + '_' + N5 + '.sorted.bam')
    elif filetype == 'sorted_bai_genome_global':
        return os.path.join(main_folder, 'bam_genome_global_files', N7 + '_' + N5 + '.sorted.bam.bai')
    elif filetype == 'results_amplicons':
        return os.path.join(main_folder, 'results', N7 + '_' + N5)  # We will append the window size later
    elif filetype == 'results_plasmid':
        return os.path.join(main_folder, 'results', N7 + '_' + N5 + '_results_plasmid.xlsx')
    elif filetype == 'results_all_amplicons':
        return os.path.join(main_folder, 'results', N7 + '_' + N5 + '_results_all_amplicons.xlsx')
    elif filetype == 'results_genomewide':
        return os.path.join(main_folder, 'results', N7 + '_' + N5 + '_results_genomewide.xlsx')
    elif filetype == 'summary_all_alignments':
        return os.path.join(main_folder, 'results', N7 + '_' + N5 + '_summary_all_alignments.xlsx')
    elif filetype == 'read_counts':
        return os.path.join(main_folder, 'results', N7 + '_' + N5 + '_read_counts.xlsx')
    ## more added
    elif filetype == 'results_pipeline2_global':
        return os.path.join(main_folder, 'results', N7 + '_' + N5 + '_results_pipeline2_global.xlsx')
    elif filetype == 'results_pipeline2_local':
        return os.path.join(main_folder, 'results', N7 + '_' + N5 + '_results_pipeline2_local.xlsx')
    elif filetype == 'sam_genome_global_single':
        return os.path.join(main_folder, 'sam_genome_global_files', N7 + '_' + N5 + '_single.sam')
    elif filetype == 'sam_genome_local_single':
        return os.path.join(main_folder, 'sam_genome_local_files', N7 + '_' + N5 + '_single.sam')
    elif filetype == 'sam_report_genome_global_single':
        return os.path.join(main_folder, 'sam_genome_global_files', N7 + '_' + N5 + '_single.sam.report.txt')
    elif filetype == 'sam_report_genome_local_single':
        return os.path.join(main_folder, 'sam_genome_local_files', N7 + '_' + N5 + '_single.sam.report.txt')
    elif filetype == 'bam_genome_global_single':
        return os.path.join(main_folder, 'bam_genome_global_files', N7 + '_' + N5 + '_single.bam')
    elif filetype == 'bam_genome_local_single':
        return os.path.join(main_folder, 'bam_genome_local_files', N7 + '_' + N5 + '_single.bam')
    elif filetype == 'sorted_bam_genome_global_single':
        return os.path.join(main_folder, 'bam_genome_global_files', N7 + '_' + N5 + '_single.sorted.bam')
    elif filetype == 'sorted_bam_genome_local_single':
        return os.path.join(main_folder, 'bam_genome_local_files', N7 + '_' + N5 + '_single.sorted.bam')
    elif filetype == 'sam_genome_global_paired':
        return os.path.join(main_folder, 'sam_genome_global_files', N7 + '_' + N5 + '_paired.sam')
    elif filetype == 'sam_genome_local_paired':
        return os.path.join(main_folder, 'sam_genome_local_files', N7 + '_' + N5 + '_paired.sam')
    elif filetype == 'sam_report_genome_global_paired':
        return os.path.join(main_folder, 'sam_genome_global_files', N7 + '_' + N5 + '_paired.sam.report.txt')
    elif filetype == 'sam_report_genome_local_paired':
        return os.path.join(main_folder, 'sam_genome_local_files', N7 + '_' + N5 + '_paired.sam.report.txt')
    elif filetype == 'bam_genome_global_paired':
        return os.path.join(main_folder, 'bam_genome_global_files', N7 + '_' + N5 + '_paired.bam')
    elif filetype == 'bam_genome_local_paired':
        return os.path.join(main_folder, 'bam_genome_local_files', N7 + '_' + N5 + '_paired.bam')
    elif filetype == 'sorted_bam_genome_global_paired':
        return os.path.join(main_folder, 'bam_genome_global_files', N7 + '_' + N5 + '_paired.sorted.bam')
    elif filetype == 'sorted_bam_genome_local_paired':
        return os.path.join(main_folder, 'bam_genome_local_files', N7 + '_' + N5 + '_paired.sorted.bam')
    elif filetype == 'filtered_and_sorted_genome_global_single':
        return os.path.join(main_folder, 'bam_genome_global_files', N7 + '_' + N5 + '_single_primary_as_mapq.sorted.bam')    
    elif filetype == 'filtered_and_sorted_genome_local_single':
        return os.path.join(main_folder, 'bam_genome_local_files', N7 + '_' + N5 + '_single_primary_as_mapq.sorted.bam')    
    elif filetype == 'genome_global_bed_single':
        return os.path.join(main_folder, 'bam_genome_global_files', N7 + '_' + N5 + '_single_global.sorted.bed')
    elif filetype == 'filtered_and_sorted_genome_global_paired':
        return os.path.join(main_folder, 'bam_genome_global_files', N7 + '_' + N5 + '_paired_primary_as_mapq.sorted.bam')  
    elif filetype == 'filtered_and_sorted_genome_local_paired':
        return os.path.join(main_folder, 'bam_genome_local_files', N7 + '_' + N5 + '_paired_primary_as_mapq.sorted.bam')  
    elif filetype == 'genome_global_bed_paired':
        return os.path.join(main_folder, 'bam_genome_global_files', N7 + '_' + N5 + '_paired_global.sorted.bed')
    elif filetype == 'guideseq_global_bed_single':
        return os.path.join(main_folder, 'bam_genome_global_files', N7 + '_' + N5 + '_single.global.sorted.bed')
    elif filetype == 'guideseq_global_bed_paired':
        return os.path.join(main_folder, 'bam_genome_global_files', N7 + '_' + N5 + '_paired.global.sorted.bed')
    elif filetype == 'guideseq_local_bed_single':
        return os.path.join(main_folder, 'bam_genome_global_files', N7 + '_' + N5 + '_single.local.sorted.bed')
    elif filetype == 'guideseq_local_bed_paired':
        return os.path.join(main_folder, 'bam_genome_global_files', N7 + '_' + N5 + '_paired.local.sorted.bed')
###########################
#
#   my homemade function to make a fastq with only correct priming events

#   the sequence to test for guideSeq is Read1
#    It trims off the primer binding sequence and downstream sequence for good reads (primer_seq_plus_downstream)
#######################

def correct_priming_guideseq(dir_sample, amplicon_info, primer_seq, primer_seq_plus_downstream):
    
    #defined the sequence as input (primer and 12 nt downstream). Normally 32-37nt sequence
    
    
    length_primer_down = len(primer_seq_plus_downstream)
    length_primer = len(primer_seq)
    
    '''
    files_out = list()
    n_file = 0
    files_out_dict = dict()

    for file_selected in files_out:
        files_out_dict[os.path.basename(file_selected)] = n_file
        n_file += 1
    '''
    
    #define the input files and output files
    N7 = amplicon_info['index_I1']
    N5 = amplicon_info['index_I2']
    name = amplicon_info['name']
    
    #input files
    r1_fastq = create_filename(dir_sample, N7, N5, 'R1fastqgz')
    r2_fastq = create_filename(dir_sample, N7, N5, 'R2fastqgz')
    
    #output files not gz compressed
    file_R1_corrprime = create_filename(dir_sample, N7, N5, 'R1fastq_CorrPrime')
    file_R2_corrprime = create_filename(dir_sample, N7, N5, 'R2fastq_CorrPrime')

    #mismatched files
    file_out_r1_not_corrprime = os.path.join(dir_sample, N7 + '_' + N5, 'fastq_files', 'mis_prime_R1.fastq')
    file_out_r2_not_corrprime = os.path.join(dir_sample, N7 + '_' + N5, 'fastq_files', 'mis_prime_R2.fastq')
    

    
    #open/create all these files
    ref_file_R1_corrprime = open(file_R1_corrprime, "w")
    ref_file_R2_corrprime = open(file_R2_corrprime, "w")

    ref_file_out_r1_not_corrprime = open(file_out_r1_not_corrprime, "w")
    ref_file_out_r2_not_corrprime = open(file_out_r2_not_corrprime, "w")
    
    
    # We open r1,r2 files and distribute reads
    with open_fastq_or_gz(r1_fastq) as r1_file, open_fastq_or_gz(r2_fastq) as r2_file:
        # Add counters for all reads

        total_reads = 0

        guideseq_primer_read = 0

        guideseq_primer_read_plusdownstream = 0
    
 
        r1_r2 = itertools.izip(r1_file, r2_file)

        for header_r1, header_r2 in r1_r2:
            total_reads += 1
            
            seq_r1, seq_r2 = r1_r2.next()

            r1_r2.next()

            qual_r1, qual_r2 = r1_r2.next()
            seq_r1, seq_r2 = seq_r1.rstrip(), seq_r2.rstrip()
            qual_r1, qual_r2 = qual_r1.rstrip(), qual_r2.rstrip()

            #We mask with N any bases with scores below or equal to , (11, default in mask)
            
            seq_r1_primer_size = mask(seq_r1[:length_primer], qual_r1[:length_primer])
            
            seq_r1_primer_and_downstream = mask(seq_r1[:length_primer_down], qual_r1[:length_primer_down])
            
            #make seq1_r1 and qual_r1 without the primer seq
            seq_r1_noprimer = seq_r1[length_primer_down:]
            qual_r1_noprimer = qual_r1[length_primer_down:]

            #this is modified for miniseq UMI->Index2
            #seq_r2_use = mask(seq_r2, qual_r1)

            # change to 1 for reads with perfect match or match within hamming distance decided above
            is_good_read = 0
            
            #tests if there is a match for the primer and the nucletides downstream
            if seq_r1_primer_and_downstream == primer_seq_plus_downstream:
                # perfect match case
                is_good_read = 1
            #tests if there is a match for the primer (could have or not half downstream)
            if seq_r1_primer_size == primer_seq:
                guideseq_primer_read += 1
                        
            #Print good reads 
            if is_good_read:
                guideseq_primer_read_plusdownstream += 1
                
                # We test whether the read has on of the combination of indices from our experiment list
                # If not save in a separate file
                r1f = ref_file_R1_corrprime
                r2f = ref_file_R2_corrprime

                print("\n".join([header_r1.rstrip(), seq_r1_noprimer.rstrip(), "+", qual_r1_noprimer.rstrip()]), file=r1f)
                print("\n".join([header_r2.rstrip(), seq_r2.rstrip(), "+", qual_r2.rstrip()]), file=r2f)
                    
        
            else:
                # We print reads with mispriming
                print("\n".join([header_r1.rstrip(), seq_r1.rstrip(), "+", qual_r1.rstrip()]), file=ref_file_out_r1_not_corrprime)
                print("\n".join([header_r2.rstrip(), seq_r2.rstrip(), "+", qual_r2.rstrip()]), file=ref_file_out_r2_not_corrprime)
                
        
    results_file = create_filename(dir_sample, N7, N5, 'localpriming_data')
    
    
    results_df = pd.DataFrame({'sample_name': [name],
                               'i7': [N7],
                               'i5': [N5],
                               'total_reads': [total_reads],
                               'reads_with_good_priming': [guideseq_primer_read_plusdownstream],
                               'reads_with_guideseq_primer_misprimed': [guideseq_primer_read - guideseq_primer_read_plusdownstream],
                               'reads_with_indexing_primer_mispriming': [total_reads - guideseq_primer_read], #indexing primer i
                               },
                                  columns=['sample_name',
                                           'i7',
                                           'i5',
                                           'total_reads',
                                           'reads_with_good_priming',
                                           'reads_with_guideseq_primer_misprimed',
                                           'reads_with_indexing_primer_mispriming'])
    
    results_df.to_excel(results_file)
    return results_df

    
############################
#
# #This was changed as it requres a dfiferent adapter sequence
#
#Remove adapters in fastq files. The idea of this is that if the sequence runs beyond the length of the acutal genomice sequence into the sequencing
#       primers on the other other side, it will then be trimmed down so you dont try and align adapter sequences. 
# Input: directory to be analyzed (fastq files)
#        dir_sample, name of the directory for the whole run, typically with the name of a miseq run
#        amplicon_info, slice of sample_info.csv for the sample being processed to get the indexes required for saving names
#        process_AMP_seq_run, set to 1 to trim in read2 the same adapter as in GUIDE-Seq
     # it is very important! to pay attention to if youre using the nextera or trueseq adapters as the sequence to trim will be different on the i7 side.
#

#
# ##########################
def trim_fastq(dir_sample, amplicon_info, process_AMP_seq_run=0):

    # UDiTaS adapters
    Nv2F = 'TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG' #for the i5 side end of read2 for tn5
    SBS12nextera = 'GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG'  #this is for the i7 side for nextera end of read1 for tnt5
  #  SBS12 = 'GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT'  #this is for i7 side for trueseq primers **need this one for the pytest**

    
    if process_AMP_seq_run == 1:
        #the primer needed for LAM
        i2_adapter = 'A'
        #i2_adapter = 'ACACTCTTTCCCTACACGACGCTCTTCCGATCT'
    else:
        i2_adapter = Nv2F

    # We first check if the experiment had any guides
    N7 = amplicon_info['index_I1']
    N5 = amplicon_info['index_I2']
    
    file_R1 = create_filename(dir_sample, N7, N5, 'R1fastq_CorrPrime')
    file_R2 = create_filename(dir_sample, N7, N5, 'R2fastq_CorrPrime')

    file_cutadapt_R1 = create_filename(dir_sample, N7, N5, 'R1trimmed')
    file_cutadapt_R2 = create_filename(dir_sample, N7, N5, 'R2trimmed')
    file_cutadapt_report = create_filename(dir_sample, N7, N5, 'trimmed_report')

    if not os.path.exists(os.path.dirname(file_cutadapt_R1)):
        os.mkdir(os.path.dirname(file_cutadapt_R1))

    # remove adapters with cutadapt
    #original uditas peramiter had an error -e 0.33 (but was cutting of random stuff too much)
    cutadapt_command = ['cutadapt',
                        '-m', '10',
                        '-e', '0.20',
                        '-a', reverse_complement(SBS12nextera),
                        '-A', reverse_complement(i2_adapter),
                        '-o', file_cutadapt_R1, '-p', file_cutadapt_R2,
                        file_R1, file_R2]

    handle_cutadapt_report = open(file_cutadapt_report, 'wb')
    subprocess.call(cutadapt_command, stdout=handle_cutadapt_report)
    handle_cutadapt_report.close()


#################################################################################
# Function to write reference amplicons with various structural rearrangements
#################################################################################
def write_amplicon(dir_sample, amplicon_info, amplicon_list):
    N7 = amplicon_info['index_I1']
    N5 = amplicon_info['index_I2']
    exp_dir = create_filename(dir_sample, N7, N5, 'mainfolder')
    amplicon_folder = os.path.join(exp_dir, 'amplicons')
    if not os.path.exists(amplicon_folder):
        os.mkdir(amplicon_folder)

    filename = os.path.join(exp_dir, amplicon_folder, 'amplicons.fa')
    file_handle = open(filename, "w")

    for amps in amplicon_list:
        seq1 = Seq(amps[1], IUPAC.unambiguous_dna)
        record1 = SeqRecord(seq1, amps[0], description='')
        SeqIO.write(record1, file_handle, 'fasta')

    file_handle.close()
    # Create index file
    initial_dir = os.getcwd()
    os.chdir(amplicon_folder)
    index_err_file = os.path.join(amplicon_folder, 'index.err')
    index_out_file = os.path.join(amplicon_folder, 'index.out')

    index_err_fh = open(index_err_file, 'wb')
    index_out_fh = open(index_out_file, 'wb')
    subprocess.call(['bowtie2-build',
                     filename, 'amplicons'], stderr=index_err_fh, stdout=index_out_fh)
    os.chdir(initial_dir)
    index_err_fh.close()
    index_out_fh.close()


####make sure there is no redundant sequences between your targeting vector and the delivery plasmid####

def align_plasmid_local(dir_sample, amplicon_info, ncpu=4):

    # We first check if the experiment had any guides
    N7 = amplicon_info['index_I1']
    N5 = amplicon_info['index_I2']
    # exp_dir = create_filename(dir_sample, N7, N5, 'mainfolder')

    file_cutadapt_R1 = create_filename(dir_sample, N7, N5, 'R1trimmed')
    file_cutadapt_R2 = create_filename(dir_sample, N7, N5, 'R2trimmed')

    file_sam_plasmid_local = create_filename(dir_sample, N7, N5, 'sam_plasmid_local')
    file_sam_report_plasmid_local = create_filename(dir_sample, N7, N5, 'sam_report_plasmid_local')

    if not os.path.exists(os.path.dirname(file_sam_plasmid_local)):
        os.mkdir(os.path.dirname(file_sam_plasmid_local))

    file_bam_plasmid_local = create_filename(dir_sample, N7, N5, 'bam_plasmid_local')
    file_sorted_bam_plasmid_local = create_filename(dir_sample, N7, N5, 'sorted_bam_plasmid_local')
    # file_sorted_bai_genome_local = create_filename(dir_sample, N7, N5, 'sorted_bai_genome_local')

    if not os.path.exists(os.path.dirname(file_bam_plasmid_local)):
        os.mkdir(os.path.dirname(file_bam_plasmid_local))

    # local alignment to the genome with bowtie2
    initial_dir = os.getcwd()

    folder_amplicons = create_filename(dir_sample, N7, N5, 'amplicons')

    os.chdir(folder_amplicons)

    bowtie2_command = ['bowtie2', '--local', '-p', str(ncpu),
                       '-X', '5000', '-k', '2', '-x', 'plasmid',
                             '-1', file_cutadapt_R1, '-2', file_cutadapt_R2,
                             '-S', file_sam_plasmid_local]

    handle_sam_report_genome_local = open(file_sam_report_plasmid_local, 'wb')

    subprocess.call(bowtie2_command, stderr=handle_sam_report_genome_local)

    handle_sam_report_genome_local.close()

    # convert sam to bam
    sam_to_bam_plasmid_local_command = ['samtools', 'view', '-Sb', file_sam_plasmid_local]

    handle_file_bam_plasmid_local = open(file_bam_plasmid_local, 'wb')

    subprocess.call(sam_to_bam_plasmid_local_command, stdout=handle_file_bam_plasmid_local)

    # sort bam files
    sort_bam_plasmid_local_command = ['samtools', 'sort', file_bam_plasmid_local, '-o', file_sorted_bam_plasmid_local]

    subprocess.call(sort_bam_plasmid_local_command)

    # Create bam index files
    create_bam_plasmid_local_index_command = ['samtools', 'index', file_sorted_bam_plasmid_local]
    subprocess.call(create_bam_plasmid_local_index_command)

    # Clean up
    os.remove(file_sam_plasmid_local)
    os.remove(file_bam_plasmid_local)

    os.chdir(initial_dir)

def extract_unmapped_reads_plasmid(dir_sample, amplicon_info):

    N7 = amplicon_info['index_I1']
    N5 = amplicon_info['index_I2']

    file_sorted_bam_plasmid_local = create_filename(dir_sample, N7, N5, 'sorted_bam_plasmid_local')

    file_unmapped_bam_plasmid = create_filename(dir_sample, N7, N5, 'unmapped_bam_plasmid_local')

    file_qsorted_unmapped_bam_plasmid = create_filename(dir_sample, N7, N5, 'qsorted_unmapped_bam_plasmid_local')

    file_R1_unmapped = create_filename(dir_sample, N7, N5, 'unmapped_plasmid_R1fastq')
    file_R2_unmapped = create_filename(dir_sample, N7, N5, 'unmapped_plasmid_R2fastq')

    if not os.path.exists(os.path.dirname(file_R1_unmapped)):
        os.mkdir(os.path.dirname(file_R1_unmapped))

    extract_unmapped_bam_command = ['samtools', 'view', '-b', '-f', '0x4', file_sorted_bam_plasmid_local, '-o',
                                    file_unmapped_bam_plasmid]

    subprocess.call(extract_unmapped_bam_command)

    qsort_unmapped_bam_command = ['samtools', 'sort', '-n', file_unmapped_bam_plasmid, '-o',
                                  file_qsorted_unmapped_bam_plasmid]

    subprocess.call(qsort_unmapped_bam_command)

    bamtofastq_command = ['bedtools', 'bamtofastq', '-i', file_qsorted_unmapped_bam_plasmid,
                          '-fq', file_R1_unmapped, '-fq2', file_R2_unmapped]

    file_err = file_R1_unmapped[:-9] + '_err.txt'
    handle_file_err = open(file_err, 'wb')

    subprocess.call(bamtofastq_command, stderr=handle_file_err)

    for fo in [file_R1_unmapped, file_R2_unmapped]:
        with open(fo) as f_in, gzip.open(fo + '.gz', 'wb') as f_out:
            f_out.writelines(f_in)
        os.remove(fo)

#analyze plasmid alignments

def analyze_alignments_plasmid(dir_sample, amplicon_info, min_MAPQ, file_genome_2bit, do_plasmid):
    N7 = amplicon_info['index_I1']
    N5 = amplicon_info['index_I2']
        
    exp_dir = create_filename(dir_sample, N7, N5, 'mainfolder')

    file_UMI = create_filename(dir_sample, N7, N5, 'umifastqgz')
    UMI_dict = create_barcode_dict(file_UMI)
    
    results_folder = os.path.join(exp_dir, 'results')
    if not os.path.exists(results_folder):
        os.mkdir(results_folder)

    results_file = create_filename(dir_sample, N7, N5, 'results_plasmid')

    if do_plasmid:
        file_sorted_bam_plasmid_local = create_filename(dir_sample, N7, N5, 'sorted_bam_plasmid_local')

        bam_in_alignment_file = pysam.AlignmentFile(file_sorted_bam_plasmid_local, 'rb')
        bam_in = bam_in_alignment_file.fetch()

        genome = twobitreader.TwoBitFile(file_genome_2bit)  # Load genome. Used for getting the sequences
        
        length_to_test = 15  # We check this number of bases after the primer
        uditas_primer_length = amplicon_info['end'] - amplicon_info['start']
        
        if amplicon_info['strand'] == '+':  # This is the UDiTaS oligo strand
            #I had to add int() command to make this work for some reason
            seq_after_uditas_primer = genome[amplicon_info['chr']][int(amplicon_info['end']):int((amplicon_info['end'] + length_to_test))]
            
        elif amplicon_info['strand'] == '-':
            seq_after_uditas_primer = reverse_complement(genome[amplicon_info['chr']][int((amplicon_info['start'] - length_to_test)):(int(amplicon_info['start']))])
        n_max_mismatches = 2  # We allow this number of mismatches between the read and the sequence after the primer

        names_list_plasmid_genome = []
        UMI_list_plasmid_genome = []
        names_list_plasmid_only = []
        UMI_list_plasmid_only = []
        
        for read in bam_in:
            if read.mapping_quality >= min_MAPQ and not read.is_unmapped and not read.is_secondary:
                if read.is_read2:  # R2 is the UDiTaS primer
                    if read.is_reverse:
                        seq_test = reverse_complement(read.query_sequence)[int(uditas_primer_length):int((uditas_primer_length + length_to_test))]
                    else:
                        seq_test = read.query_sequence[int(uditas_primer_length): int(uditas_primer_length + length_to_test)]
                    # Sometimes, after cutadapt we have a read shorter than uditas_primer_length + length_to_test
                    # We skip those directly without calculating hamm_dist, which doesn't make sense
                    if (len(seq_test) == len(seq_after_uditas_primer.upper()) and
                        hamm_dist(seq_test, seq_after_uditas_primer.upper()) <= n_max_mismatches):
                        # Reads for which the R2 has genomic sequence after the UDiTaS primer
                        UMI_list_plasmid_genome.append(UMI_dict[read.query_name][0])
                        names_list_plasmid_genome.append(read.query_name)
                    else: # We put those short reads into the plasmid only bucket
                        UMI_list_plasmid_only.append(UMI_dict[read.query_name][0])
                        names_list_plasmid_only.append(read.query_name)

        total_reads_plasmid_genome = len(set(names_list_plasmid_genome))
        total_reads_collapsed_plasmid_genome = len(set(UMI_list_plasmid_genome))
        total_reads_plasmid_only = len(set(names_list_plasmid_only))
        total_reads_collapsed_plasmid_only = len(set(UMI_list_plasmid_only))

        results_df = pd.DataFrame({'target_plus_plasmid_total_reads': [total_reads_plasmid_genome],
                                   'target_plus_plasmid_total_reads_collapsed': [total_reads_collapsed_plasmid_genome],
                                   'plasmid_only_total_reads': [total_reads_plasmid_only],
                                   'plasmid_only_total_reads_collapsed': [total_reads_collapsed_plasmid_only]
                                   },
                                  columns=['target_plus_plasmid_total_reads',
                                           'target_plus_plasmid_total_reads_collapsed',
                                           'plasmid_only_total_reads',
                                           'plasmid_only_total_reads_collapsed'])
    else:
        results_df = pd.DataFrame(index=np.arange(1),
                                  columns=['target_plus_plasmid_total_reads',
                                           'target_plus_plasmid_total_reads_collapsed',
                                           'plasmid_only_total_reads',
                                           'plasmid_only_total_reads_collapsed'])

    results_df.to_excel(results_file)

    return results_df

#################################################################################
# Function to create barcode dict
#################################################################################
def create_barcode_dict(filename):
    barcode_file = open_fastq_or_gz(filename)

    barcode_dict = dict()

    barcode_reads = itertools.izip(barcode_file)

    for header_barcode in barcode_reads:
        seq_barcode = barcode_reads.next()
        barcode_reads.next()
        qual_barcode = barcode_reads.next()
        barcode_dict[header_barcode[0].split()[0][1:]] = [seq_barcode[0].rstrip(), qual_barcode[0].rstrip()]

    return barcode_dict

############################
#
# Aligns reads globally to amplicon. "end-to-end" as the default function in bowtie2
# Input: directory to be analyzed
#        amplicon_info, slice of sample_info.csv for the sample being processed
#        file_genome_2bit, 2bit file with the reference genome being used
#
# ##########################
def align_amplicon(dir_sample, amplicon_info, check_plasmid_insertions, ncpu=4):

    # We first check if the experiment had any guides
    N7 = amplicon_info['index_I1']
    N5 = amplicon_info['index_I2']

    has_plasmid = type(amplicon_info['plasmid_sequence']) is str or type(amplicon_info['plasmid_sequence']) is unicode

    if check_plasmid_insertions == 1 and has_plasmid:
        file_R1 = create_filename(dir_sample, N7, N5, 'unmapped_plasmid_R1fastqgz')
        file_R2 = create_filename(dir_sample, N7, N5, 'unmapped_plasmid_R2fastqgz')
    else:
        file_R1 = create_filename(dir_sample, N7, N5, 'R1trimmed')
        file_R2 = create_filename(dir_sample, N7, N5, 'R2trimmed')

    if not os.path.exists(os.path.dirname(file_R1)):
        os.mkdir(os.path.dirname(file_R1))

    file_sam_amplicons = create_filename(dir_sample, N7, N5, 'sam_amplicons')
    file_sam_report_amplicons = create_filename(dir_sample, N7, N5, 'sam_report_amplicons')

    if not os.path.exists(os.path.dirname(file_sam_amplicons)):
        os.mkdir(os.path.dirname(file_sam_amplicons))

    file_bam_amplicons = create_filename(dir_sample, N7, N5, 'bam_amplicons')
    file_sorted_bam_amplicons = create_filename(dir_sample, N7, N5, 'sorted_bam_amplicons')

    if not os.path.exists(os.path.dirname(file_bam_amplicons)):
        os.mkdir(os.path.dirname(file_bam_amplicons))

    # global alignment to the amplicons with bowtie2
    initial_dir = os.getcwd()
    folder_amplicons = create_filename(dir_sample, N7, N5, 'amplicons')

    os.chdir(folder_amplicons)
    bowtie2_command = ['bowtie2', '-p', str(ncpu), '--very-sensitive',
                       '-X', '5000', '-k', '2', '-x', 'amplicons',
                       '-1', file_R1, '-2', file_R2,
                       '-S', file_sam_amplicons]

    handle_sam_report_amplicons = open(file_sam_report_amplicons, 'wb')

    subprocess.call(bowtie2_command, stderr=handle_sam_report_amplicons)

    handle_sam_report_amplicons.close()

    # convert sam to bam
    sam_to_bam_amplicons_command = ['samtools', 'view', '-Sb', file_sam_amplicons]

    handle_file_bam_amplicons = open(file_bam_amplicons, 'wb')

    subprocess.call(sam_to_bam_amplicons_command, stdout=handle_file_bam_amplicons)

    # sort bam files
    sort_bam_amplicons_command = ['samtools', 'sort', file_bam_amplicons, '-o', file_sorted_bam_amplicons]

    subprocess.call(sort_bam_amplicons_command)

    # Clean up
    os.remove(file_sam_amplicons)
    os.remove(file_bam_amplicons)

    # Create bam index files
    create_bam_amplicons_index_command = ['samtools', 'index', file_sorted_bam_amplicons]
    subprocess.call(create_bam_amplicons_index_command)

    os.chdir(initial_dir)

################################################################################
# Function to extract all unmapped reads to the amplicons
################################################################################
def extract_unmapped_reads_amplicons(dir_sample, amplicon_info):

    N7 = amplicon_info['index_I1']
    N5 = amplicon_info['index_I2']

    file_sorted_bam_amplicons = create_filename(dir_sample, N7, N5, 'sorted_bam_amplicons')

    file_unmapped_bam_amplicons = create_filename(dir_sample, N7, N5, 'unmapped_bam_amplicons')

    file_qsorted_unmapped_bam_amplicons = create_filename(dir_sample, N7, N5, 'qsorted_unmapped_bam_amplicons')

    file_R1_unmapped = create_filename(dir_sample, N7, N5, 'unmapped_amplicons_R1fastq')
    file_R2_unmapped = create_filename(dir_sample, N7, N5, 'unmapped_amplicons_R2fastq')
    file_unmapped_report = create_filename(dir_sample, N7, N5, 'unmapped_amplicons_report')

    if not os.path.exists(os.path.dirname(file_R1_unmapped)):
        os.mkdir(os.path.dirname(file_R1_unmapped))

    extract_unmapped_bam_command = ['samtools', 'view', '-b', '-f', '0x4', file_sorted_bam_amplicons, '-o',
                                    file_unmapped_bam_amplicons]

    subprocess.call(extract_unmapped_bam_command)

    qsort_unmapped_bam_command = ['samtools', 'sort', '-n', file_unmapped_bam_amplicons, '-o',
                                  file_qsorted_unmapped_bam_amplicons]

    subprocess.call(qsort_unmapped_bam_command)

    bamtofastq_command = ['bedtools', 'bamtofastq', '-i', file_qsorted_unmapped_bam_amplicons,
                          '-fq', file_R1_unmapped, '-fq2', file_R2_unmapped]

    handle_unmapped_report = open(file_unmapped_report, 'wb')
    subprocess.call(bamtofastq_command, stderr=handle_unmapped_report)

    for fo in [file_R1_unmapped, file_R2_unmapped]:
        with open(fo) as f_in, gzip.open(fo + '.gz', 'wb') as f_out:
            f_out.writelines(f_in)
        os.remove(fo)



################################################################################
# Function to get the number of intersecting bases between two intervals
################################################################################
def get_intersection(region1_begin, region1_end, region2_begin, region2_end):
    list1 = range(int(region1_begin) + 1, int(region1_end) + 1)
    list2 = range(int(region2_begin) + 1, int(region2_end) + 1)
    return len(set(list1).intersection(list2))



# Functions from Uditas to analyze the alignments. These support analyze_alignment()

# create_segments() is used in analyze_fragment_sizes(). Analyze_fragment_sizes is used in analyze_alignments()
# These functions are unchanged from Uditas

################################################################################
# helper function to create list of fragment coordinates, useful to get size statistics
################################################################################
def create_segments(iter1, bam_in, min_MAPQ):
    segments = list()
    for read in iter1:
        if read.mapping_quality >= min_MAPQ and read.is_read2 and read.is_paired:
            if read.is_reverse:
                segment_start = read.reference_end + read.tlen  # Note: read.tlen is < 0
                segment_end = read.reference_end  # pysam is 0 based
            else:
                segment_start = read.reference_start  # pysam is 0 based
                segment_end = read.reference_start + read.tlen

            if segment_start < 0:  # We don't want to go below 0
                segment_start = 0

            if segment_end > segment_start:
                segments.append((bam_in.getrname(read.reference_id), segment_start,
                                 segment_end))

    return segments





# Functions from Uditas to analyze the alignments. These support analyze_alignment()

# get_cut_in_reference_amplicon_df() is used in znalyze_alignment()
# This had to be altered to get the cuts for REPLACE



############################
#
# Aligns reads to the whole genome using bowtie2 and global alignment. This is used to find mispriming events
#  and unmapped reads
#
# Input: directory to be analyzed
#        amplicon_info, slice of sample_info.csv for the sample being processed
#        assembly, name of the assembly to be used by bowtie2. It is convenient to place it in a folder especified by
#        the environmental variable BOWTIE2_INDEXES
#
# ##########################
def align_genome_global(dir_sample, amplicon_info, assembly, ncpu=4):

    # We first check if the experiment had any guides
    N7 = amplicon_info['index_I1']
    N5 = amplicon_info['index_I2']

    file_R1 = create_filename(dir_sample, N7, N5, 'unmapped_amplicons_R1fastqgz')
    file_R2 = create_filename(dir_sample, N7, N5, 'unmapped_amplicons_R2fastqgz')

    file_sam_genome_global = create_filename(dir_sample, N7, N5, 'sam_genome_global')
    file_sam_report_genome_global = create_filename(dir_sample, N7, N5, 'sam_report_genome_global')

    if not os.path.exists(os.path.dirname(file_sam_genome_global)):
        os.mkdir(os.path.dirname(file_sam_genome_global))

    file_bam_genome_global = create_filename(dir_sample, N7, N5, 'bam_genome_global')
    file_sorted_bam_genome_global = create_filename(dir_sample, N7, N5, 'sorted_bam_genome_global')

    if not os.path.exists(os.path.dirname(file_bam_genome_global)):
        os.mkdir(os.path.dirname(file_bam_genome_global))

    # global alignment to the genome with bowtie2
    initial_dir = os.getcwd()

    bowtie2_command = ['bowtie2', '--very-sensitive', '-p', str(ncpu),
                       '-X', '5000', '-k', '2', '-x', assembly,
                       '-1', file_R1, '-2', file_R2,
                       '-S', file_sam_genome_global]

    handle_sam_report_genome_global = open(file_sam_report_genome_global, 'wb')

    subprocess.call(bowtie2_command, stderr=handle_sam_report_genome_global)

    handle_sam_report_genome_global.close()

    # convert sam to bam
    sam_to_bam_genome_global_command = ['samtools', 'view', '-Sb', file_sam_genome_global]

    handle_file_bam_genome_global = open(file_bam_genome_global, 'wb')

    subprocess.call(sam_to_bam_genome_global_command, stdout=handle_file_bam_genome_global)

    # sort bam files
    sort_bam_genome_global_command = ['samtools', 'sort', file_bam_genome_global, '-o', file_sorted_bam_genome_global]

    subprocess.call(sort_bam_genome_global_command)

    # Clean up
    os.remove(file_sam_genome_global)
    os.remove(file_bam_genome_global)

    # Create bam index files
    create_bam_genome_global_index_command = ['samtools', 'index', file_sorted_bam_genome_global]
    subprocess.call(create_bam_genome_global_index_command)

    os.chdir(initial_dir)

################################################################################
# Function to analyze global alignments to the genome
# this genome analysis really needs to be done for the excel stuff I think
################################################################################
def analyze_alignments_genome_global(dir_sample, amplicon_info, min_MAPQ, min_AS,  file_genome_2bit):
    N7 = amplicon_info['index_I1']
    N5 = amplicon_info['index_I2']

    genome = twobitreader.TwoBitFile(file_genome_2bit)

    exp_dir = create_filename(dir_sample, N7, N5, 'mainfolder')

    file_UMI = create_filename(dir_sample, N7, N5, 'umifastqgz')
    UMI_dict = create_barcode_dict(file_UMI)

    results_folder = os.path.join(exp_dir, 'results')
    if not os.path.exists(results_folder):
        os.mkdir(results_folder)

    results_file = create_filename(dir_sample, N7, N5, 'results_genomewide')

    file_sorted_bam_genome_global = create_filename(dir_sample, N7, N5, 'sorted_bam_genome_global')

    bam_in_alignment_file = pysam.AlignmentFile(file_sorted_bam_genome_global, 'rb')

    bam_in_all = bam_in_alignment_file.fetch()

    names_list_genome = []
    UMI_list_genome = []
    readcount = 0
    for read in bam_in_all:
        readcount += 1
        if read.has_tag('AS'):
            read_AS = read.get_tag('AS')
        # We test first if the read is unmapped, otherwise read_AS would be undefined
        if not read.is_unmapped and (read.mapping_quality >= min_MAPQ
                                     and read_AS >= min_AS and not read.is_secondary):
            UMI_list_genome.append(UMI_dict[read.query_name][0])
            names_list_genome.append(read.query_name)
    print('the number of alignments that were attempted to align', readcount)
    names_list_target_only = []
    UMI_list_target_only = []

    fetch_window = 1000
    fetch_chr = amplicon_info['chr']

    if amplicon_info['strand'] == '+':  # This is the UDiTaS oligo strand
        fetch_start = amplicon_info['start']
        fetch_end = amplicon_info['end'] + fetch_window
        if fetch_end > len(genome[amplicon_info['chr']]):
            fetch_end = len(genome[amplicon_info['chr']])

    elif amplicon_info['strand'] == '-':
        fetch_start = amplicon_info['start'] - fetch_window
        fetch_end = amplicon_info['end']
        if fetch_start < 0:
            fetch_start = 0

    bam_in_target = bam_in_alignment_file.fetch(fetch_chr, fetch_start, fetch_end)

    for read in bam_in_target:
        if read.mapping_quality >= min_MAPQ and not read.is_unmapped and not read.is_secondary:
            UMI_list_target_only.append(UMI_dict[read.query_name][0])
            names_list_target_only.append(read.query_name)

    genomewide_total_reads = len(set(names_list_genome))
    genomewide_total_reads_collapsed = len(set(UMI_list_genome))
    genomewide_target_only_reads = len(set(names_list_target_only))
    genomewide_target_only_reads_collapsed = len(set(UMI_list_target_only))

    results_df = pd.DataFrame({'genomewide_total_reads': [genomewide_total_reads],
                               'genomewide_total_reads_collapsed': [genomewide_total_reads_collapsed],
                               'genomewide_target_only_reads': [genomewide_target_only_reads],
                               'genomewide_target_only_reads_collapsed': [genomewide_target_only_reads_collapsed]
                               },
                              columns=['genomewide_total_reads',
                                       'genomewide_total_reads_collapsed',
                                       'genomewide_target_only_reads',
                                       'genomewide_target_only_reads_collapsed'])

    results_df.to_excel(results_file)

    return results_df

#This is a function to count the amount of original files
# I did not modify this from uditas

################################################################################
# Fast way to count reads, but only works on unix
################################################################################
def wc_unix(filename):
    cat_out = subprocess.Popen(('zcat', filename), stdout=subprocess.PIPE)
    return int(int(subprocess.check_output(["wc", "-l"], stdin=cat_out.stdout).split()[0])/4.)
################################################################################
# Function to count reads
#used in uditas.py to count total input reads
################################################################################
def count_reads(dir_sample, amplicon_info):
    N7 = amplicon_info['index_I1']
    N5 = amplicon_info['index_I2']

    read_counts_file = create_filename(dir_sample, N7, N5, 'read_counts')
    file_cutadapt_R1 = create_filename(dir_sample, N7, N5, 'R1trimmed')

    rc = wc_unix(file_cutadapt_R1)
    df = pd.DataFrame({'read_count': [rc]})

    df.to_excel(read_counts_file)

    return rc

################################################################################
#  Function to get the percentages for the alignment of all reads mapped to ALL plasmid, amplicons and genomewide
################################################################################
def get_summary_all_alignments(dir_sample, amplicon_info, read_count, result_plasmid_df,
               result_reads_in_all_amplicons_df, results_genome_global_df):
    N7 = amplicon_info['index_I1']
    N5 = amplicon_info['index_I2']

    summary_all_alignments_file = create_filename(dir_sample, N7, N5, 'summary_all_alignments')

    summary_all_alignments = pd.concat([read_count, result_plasmid_df,
               result_reads_in_all_amplicons_df, results_genome_global_df], axis=1)

    total_reads_list = [k for k in summary_all_alignments.keys() if str(k).endswith('_total_reads')]

    summary_all_alignments['total_aligned'] = summary_all_alignments[total_reads_list].sum(axis=1)

    total_reads_collapsed_list = [k for k in summary_all_alignments.keys() if str(k).endswith('_total_reads_collapsed')]

    summary_all_alignments['total_aligned_collapsed'] = summary_all_alignments[total_reads_collapsed_list].sum(axis=1)

    summary_all_alignments['percent_aligned'] = 100 * (summary_all_alignments['total_aligned'] /
                                                       summary_all_alignments['read_count'])
    summary_all_alignments['percent_aligned_all_amplicons'] = 100 * (summary_all_alignments['all_amplicons_total_reads'] /
                                                        summary_all_alignments['total_aligned'])
    summary_all_alignments.to_excel(summary_all_alignments_file, index=False)

    return summary_all_alignments

################################################################################
# Function to summarize counts into all amplicons
################################################################################
def summarize_results(results):

    total_reads_list = [k for k in results.keys() if str(k).endswith('_total_reads')]
    total_reads_collapsed_list = [k for k in results.keys() if str(k).endswith('_total_reads_collapsed')]

    results['total_aligned_junctions'] = results[total_reads_list].sum(axis=1)

    for k in total_reads_list:
        # We add np.finfo(float).eps to prevent dividing by 0
        results[k + '_percent'] = 100 * results[k] / (results['total_aligned_junctions'] +
                                                      np.finfo(float).eps)

    results['total_aligned_junctions_collapsed'] = results[total_reads_collapsed_list].sum(axis=1)

    for k in total_reads_collapsed_list:
        results[k + '_percent'] = 100 * results[k] / (results['total_aligned_junctions_collapsed'] +
                                                      np.finfo(float).eps)

    return results

#######################################################################################
# We pivot the table using melt for easier visualization with other tools like Tableau
#######################################################################################
def melt_results(results_summary_with_experiments):

    melt_list = [k for k in results_summary_with_experiments.keys() if
                 str(k).endswith('_total_reads_collapsed_percent')]

    frozen_list = list(results_summary_with_experiments)

    for el in melt_list:
        frozen_list.remove(el)

    results_out = pd.melt(results_summary_with_experiments,
                          value_vars=melt_list,
                          id_vars=frozen_list, var_name='Type', value_name='Percent Editing')

    return results_out

#######################################################################################
# We pivot but don't use the collaped reads because the UMI information doesn't exist for LAM
#######################################################################################
def melt_results_lam(results_summary_with_experiments):

    melt_list = [k for k in results_summary_with_experiments.keys() if
                 str(k).endswith('_total_reads_percent')]

    frozen_list = list(results_summary_with_experiments)

    for el in melt_list:
        frozen_list.remove(el)

    results_out = pd.melt(results_summary_with_experiments,
                          value_vars=melt_list,
                          id_vars=frozen_list, var_name='Type', value_name='Percent Editing')

    return results_out


############################
#
# Align reads that did not map to the amplicons of interest
#
# Input: directory to be analyzed
#        amplicon_info, slice of sample_info.csv for the sample being processed
#        assembly, name of the assembly to be used by bowtie2. It is convenient to place it in a folder especified by
#        the environmental variable BOWTIE2_INDEXES. Need to re-index with the targeting vector and plasmid if using
#       lam_or_tn5: needs to be 'lam' or 'tn5'
#       keep_sam: does not delete the sam file if '1'.       
#
# ##########################

def align_afterbreak_end_to_end_genome_global(dir_sample, amplicon_info, assembly, lam_or_tn5='tn5', ncpu=4, keep_sam=0):

    # We first check if the experiment had any guides
    N7 = amplicon_info['index_I1']
    N5 = amplicon_info['index_I2']
    
    #input
    # Read1 for Tn5 Read2 for LAM
    if lam_or_tn5 == 'tn5':
        file_Read = create_filename(dir_sample, N7, N5, 'R1trimmed2break')
        print('tn5 analysis')
    elif lam_or_tn5 == 'lam':
        file_Read = create_filename(dir_sample, N7, N5, 'R2trimmed2break')
        print('lam analysis')


    file_sam_genome_global = create_filename(dir_sample, N7, N5, 'break_trimmed_sam_genome_global')
    file_sam_report_genome_global = create_filename(dir_sample, N7, N5, 'break_trimmed_sam_report_genome_global')

    if not os.path.exists(os.path.dirname(file_sam_genome_global)):
        os.mkdir(os.path.dirname(file_sam_genome_global))

    file_bam_genome_global = create_filename(dir_sample, N7, N5, 'break_trimmed_bam_genome_global')
    file_sorted_bam_genome_global = create_filename(dir_sample, N7, N5, 'break_trimmed_sorted_bam_genome_global')

    if not os.path.exists(os.path.dirname(file_bam_genome_global)):
        os.mkdir(os.path.dirname(file_bam_genome_global))

    # global alignment to the genome with bowtie2
    initial_dir = os.getcwd()

    bowtie2_command = ['bowtie2', '--very-sensitive', '-p', str(ncpu),
                       '-X', '5000', '-k', '2', '-x', assembly,
                       '-U', file_Read, '-S', file_sam_genome_global]

    handle_sam_report_genome_global = open(file_sam_report_genome_global, 'wb')

    subprocess.call(bowtie2_command, stderr=handle_sam_report_genome_global)

    handle_sam_report_genome_global.close()

    # convert sam to bam
    sam_to_bam_genome_global_command = ['samtools', 'view', '-Sb', file_sam_genome_global]

    handle_file_bam_genome_global = open(file_bam_genome_global, 'wb')

    subprocess.call(sam_to_bam_genome_global_command, stdout=handle_file_bam_genome_global)

    # sort bam files
    sort_bam_genome_global_command = ['samtools', 'sort', file_bam_genome_global, '-o', file_sorted_bam_genome_global]

    subprocess.call(sort_bam_genome_global_command)

    # Clean up
    if keep_sam == 0:
        os.remove(file_sam_genome_global)
        print('sam file deleted')

    
    os.remove(file_bam_genome_global)

    # Create bam index files
    create_bam_genome_global_index_command = ['samtools', 'index', file_sorted_bam_genome_global]
    subprocess.call(create_bam_genome_global_index_command)

    os.chdir(initial_dir)


############################
#       If you put in unmapped = 1 then it uses the unmapped to amplicons read, otherwise it just maps from trimmed fastq
# Aligns reads to the whole genome using bowtie2 and local alignment. This is needed to find translocations using
#   split reads and a program like socrates
#   --local aligment means that the sequence can be soft clipped at ends and is used in translocation programs that then align the clipped ends.
#
# Input: directory to be analyzed
#        amplicon_info, slice of sample_info.csv for the sample being processed
#        assembly, name of the assembly to be used by bowtie2. It is convenient to place it in a folder especified by
#        the environmental variable BOWTIE2_INDEXES (make sure it is re-indexed with targeting vector and plasmid)
#       keep_sam: this will not delete the sam file after alignment
#       
#
# ##########################
def align_afterbreaks_genome_local(dir_sample, unmapped_amplicons, amplicon_info, assembly, lam_or_tn5='tn5', ncpu=4, keep_sam=0):

    # We first check if the experiment had any guides
    N7 = amplicon_info['index_I1']
    N5 = amplicon_info['index_I2']
    

    
    #input
    # Read1 for Tn5 Read2 for LAM
    if lam_or_tn5 == 'tn5':
        file_Read = create_filename(dir_sample, N7, N5, 'R1trimmed2break')
        print('tn5 analysis')
    elif lam_or_tn5 == 'lam': 
        file_Read = create_filename(dir_sample, N7, N5, 'R2trimmed2break')
        print('lam analysis')
    
    file_sam_genome_local = create_filename(dir_sample, N7, N5, 'break_trimmed_sam_genome_local')
    file_sam_report_genome_local = create_filename(dir_sample, N7, N5, 'break_trimmed_sam_report_genome_local')

    if not os.path.exists(os.path.dirname(file_sam_genome_local)):
        os.mkdir(os.path.dirname(file_sam_genome_local))

    file_bam_genome_local = create_filename(dir_sample, N7, N5, 'break_trimmed_bam_genome_local')
    file_sorted_bam_genome_local = create_filename(dir_sample, N7, N5, 'break_trimmed_sorted_bam_genome_local')
    # file_sorted_bai_genome_local = create_filename(dir_sample, N7, N5, 'sorted_bai_genome_local')

    if not os.path.exists(os.path.dirname(file_bam_genome_local)):
        os.mkdir(os.path.dirname(file_bam_genome_local))

    # local alignment to the genome with bowtie2
    initial_dir = os.getcwd()

    bowtie2_command = ['bowtie2', '--local', '-p', str(ncpu),
                       '-X', '5000', '-k', '2', '-x', assembly,
                             '-U', file_Read, '-S', file_sam_genome_local]

    handle_sam_report_genome_local = open(file_sam_report_genome_local, 'wb')

    subprocess.call(bowtie2_command, stderr=handle_sam_report_genome_local)

    handle_sam_report_genome_local.close()

    # convert sam to bam
    sam_to_bam_genome_local_command = ['samtools', 'view', '-Sb', file_sam_genome_local]

    handle_file_bam_genome_local = open(file_bam_genome_local, 'wb')

    subprocess.call(sam_to_bam_genome_local_command, stdout=handle_file_bam_genome_local)

    # sort bam files
    sort_bam_genome_local_command = ['samtools', 'sort', file_bam_genome_local, '-o', file_sorted_bam_genome_local]

    subprocess.call(sort_bam_genome_local_command)

    # Clean up
    if keep_sam == 0:
        os.remove(file_sam_genome_local)
        print('sam file deleted')
    os.remove(file_bam_genome_local)

    # Create bam index files
    create_bam_genome_local_index_command = ['samtools', 'index', file_sorted_bam_genome_local]
    subprocess.call(create_bam_genome_local_index_command)

    os.chdir(initial_dir)




#############################
#
#  This is a function to sort the alignemnt files. It takes primary alignments with good AS scores and MapQ and makes a bam file
#
#  input: directory of the file.
# global_or_local is if it is sorting o ut the bam file made by the global or local alignment
# tn5_or_lam is just for organiznig het output plot. Lam has no UMI collapsing so the table just removes the umi info

# min_MAPq is the minimum map quality in the alignment files
# min_AS this is the minimal alignment socre. Not sure what are good numbers but Uditas said -180 for global alignments
#
############################

def final_trimmed_bam_filtered_mapq_AS_primary(dir_sample, global_or_local, amplicon_info, tn5_or_lam = 'tn5',
                                               min_MAPQ = 25, min_AS=-180):
    #function beginning
    # define samples
    N7 = amplicon_info['index_I1']
    N5 = amplicon_info['index_I2']

    file_UMI = create_filename(dir_sample, N7, N5, 'umifastqgz')
    UMI_dict = create_barcode_dict(file_UMI)

    
    if global_or_local == 'global':
        bam_file = create_filename(dir_sample, N7, N5, 'break_trimmed_sorted_bam_genome_global')
        final_trimmed_bam_filtered_file = create_filename(dir_sample, N7, N5, 'break_trimmed_filtered_and_sorted_genome_global')
    elif global_or_local == 'local':
        bam_file = create_filename(dir_sample, N7, N5, 'break_trimmed_sorted_bam_genome_local')
        final_trimmed_bam_filtered_file = create_filename(dir_sample, N7, N5, 'break_trimmed_filtered_and_sorted_bam_genome_local')
    
    all_reads_unfiltered = []
    all_reads_unfiltered_umis = []
    all_reads_pass_filter = []
    all_reads_pass_filter_umis = []
    
    bam_alignment_file = pysam.AlignmentFile(bam_file, 'rb')
    bam_in = bam_alignment_file.fetch()
    
    filtered_reads = pysam.AlignmentFile(final_trimmed_bam_filtered_file, "wb", template=bam_alignment_file)
    
    #test each read
    for read in bam_in:
        all_reads_unfiltered.append(read.query_name)
        all_reads_unfiltered_umis.append(UMI_dict[read.query_name][0])
        if read.has_tag('AS'):
            read_AS = read.get_tag('AS')
        #check read quality is good
        if read.mapping_quality >= min_MAPQ and read_AS >= min_AS and  not read.is_secondary:
            all_reads_pass_filter.append(read.query_name)
            all_reads_pass_filter_umis.append(UMI_dict[read.query_name][0])
            filtered_reads.write(read)
    filtered_reads.close()
    bam_alignment_file.close()


    count_unfiltered = len(all_reads_unfiltered)
    count_unfiltered_umis = len(set(all_reads_unfiltered_umis))

    count_filtered = len(all_reads_pass_filter)
    count_filtered_umis = len(set(all_reads_pass_filter_umis))
    
    results_df = pd.DataFrame({'all_alignments_count': [count_unfiltered],
                               'all_alignments_count_collapsed': [count_unfiltered_umis],
                               'filtered_alignments_count': [count_filtered],
                               'filtered_alignments_count_collapsed': [count_filtered_umis],
                              },
                                  columns = ['all_alignments_count',
                                             'all_alignments_count_collapsed', 
                                             'filtered_alignments_count',
                                             'filtered_alignments_count_collapsed'])
    
    if tn5_or_lam == 'lam':
        results_df = results_df[['all_alignments_count','filtered_alignments_count']]
    
    return results_df

############################
#
# #This was changed for guideSeq they way we did it with chu
# #I already removed the 5' primer sequence on Read1 during good priming filtering
#  This is to remove the sequence for hsort reads that overlaps
#    
# Input: directory to be analyzed (fastq files)
#        dir_sample, name of the directory for the whole run, typically with the name of a miseq run
#        amplicon_info, slice of sample_info.csv for the sample being processed to get the indexes required for saving names
#        process_AMP_seq_run, set to 1 to trim in read2 the same adapter as in GUIDE-Seq



#
# ##########################


def trim_guideseq(dir_sample, amplicon_info):

    # Read1 This is the sequence for the gene specific primer binding
    R1primer = 'ATACCGTTATTAACATATGACAACTCAATTAAAC'  #for the i5 side end of read2 for tn5
    direction3 = 'GTTTAATTGAGTTGTCATATGTTAATAACGGTAT'  #this is for the i7 side for nextera end of read1 for tnt5
    # reverse trimming which is the universal Tn5 priming
    nexteraR1 ='GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG'
    nexteraR1_rev = reverse_complement(nexteraR1)
    
    direction = amplicon_info['Direction']
    if direction == 3:
        R1primer = direction3
    elif direction == 5:
        R1primer = direction5


    # We first check if the experiment had any guides
    N7 = amplicon_info['index_I1']
    N5 = amplicon_info['index_I2']
    
    file_R1 = create_filename(dir_sample, N7, N5, 'R1fastq_CorrPrime')
    file_R2 = create_filename(dir_sample, N7, N5, 'R2fastq_CorrPrime')

    file_cutadapt_R1 = create_filename(dir_sample, N7, N5, 'R1trimmed')
    file_cutadapt_R2 = create_filename(dir    if direction == 3:
        R1primer = direction3_sample, N7, N5, 'R2trimmed')
    file_cutadapt_report = create_filename(dir_sample, N7, N5, 'trimmed_report')
    
    
    if not os.path.exists(os.path.dirname(file_cutadapt_R1)):
        os.mkdir(os.path.dirname(file_cutadapt_R1))

    # remove adapters with cutadapt
    #original uditas peramiter had an error -e 0.33 (but was cutting of random stuff too much)
    # -a is hte 3' adapter for Read1
    # -A is 3' adapter for Read2
    # -m minium length
    cutadapt_command = ['cutadapt',
                        '-m', '50',
                        '-e', '0.20',
                        '-a', nexteraR1_rev,
                        '-A', reverse_complement(R1primer),
                        '-o', file_cutadapt_R1, '-p', file_cutadapt_R2,
                        file_R1, file_R2]

    handle_cutadapt_report = open(file_cutadapt_report, 'wb')
    subprocess.call(cutadapt_command, stdout=handle_cutadapt_report)
    handle_cutadapt_report.close()




############################
#
# Align reads that did not map to the amplicons of interest
#
# Input: directory to be analyzed
#        amplicon_info, slice of sample_info.csv for the sample being processed
#        assembly, name of the assembly to be used by bowtie2. It is convenient to place it in a folder especified by
#        the environmental variable BOWTIE2_INDEXES. Need to re-index with the targeting vector and plasmid if using
#       lam_or_tn5: needs to be 'lam' or 'tn5'
#       keep_sam: does not delete the sam file if '1'.  
#      single_or_paired means align single en dor paired end read     
#
# ##########################


def align_guideseq_end_to_end_genome_global(dir_sample, amplicon_info, assembly, single_or_paired = 'single', ncpu=12, keep_sam=0):

    # We first check if the experiment had any guides
    N7 = amplicon_info['index_I1']
    N5 = amplicon_info['index_I2']
    
    #input
    # Read1 for guideSeq they way Chu did it
    file_R1 = create_filename(dir_sample, N7, N5, 'R1trimmed')
    file_R2 = create_filename(dir_sample, N7, N5, 'R2trimmed')
    
    if single_or_paired == 'single':
        file_sam_genome_global = create_filename(dir_sample, N7, N5, 'sam_genome_global_single')
        file_sam_report_genome_global = create_filename(dir_sample, N7, N5, 'sam_report_genome_global_single')
    
        if not os.path.exists(os.path.dirname(file_sam_genome_global)):
            os.mkdir(os.path.dirname(file_sam_genome_global))
    
        file_bam_genome_global = create_filename(dir_sample, N7, N5, 'bam_genome_global_single')
        file_sorted_bam_genome_global = create_filename(dir_sample, N7, N5, 'sorted_bam_genome_global_single')

        if not os.path.exists(os.path.dirname(file_bam_genome_global)):
            os.mkdir(os.path.dirname(file_bam_genome_global))
    
    elif single_or_paired == 'paired':
        file_sam_genome_global = create_filename(dir_sample, N7, N5, 'sam_genome_global_paired')
        file_sam_report_genome_global = create_filename(dir_sample, N7, N5, 'sam_report_genome_global_paired')
    
        if not os.path.exists(os.path.dirname(file_sam_genome_global)):
            os.mkdir(os.path.dirname(file_sam_genome_global))
    
        file_bam_genome_global = create_filename(dir_sample, N7, N5, 'bam_genome_global_paired')
        file_sorted_bam_genome_global = create_filename(dir_sample, N7, N5, 'sorted_bam_genome_global_paired')

        if not os.path.exists(os.path.dirname(file_bam_genome_global)):
            os.mkdir(os.path.dirname(file_bam_genome_global))
    
        
    # global alignment to the genome with bowtie2
    initial_dir = os.getcwd()
    
    #-X length for paired end fragment -k number of distict matches
    if single_or_paired == 'single':
        
        bowtie2_command = ['bowtie2', '--very-sensitive', '-p', str(ncpu),
                           '-k', '4', '-x', assembly,
                           '-U', file_R1, '-S', file_sam_genome_global]
    elif single_or_paired == 'paired':
        bowtie2_command = ['bowtie2', '--very-sensitive', '-p', str(ncpu),
                            '-X', '5000', '-k', '4', '-x', assembly,
                            '-1', file_R1, '-2', file_R2,
                            '-S', file_sam_genome_global]

    handle_sam_report_genome_global = open(file_sam_report_genome_global, 'wb')

    subprocess.call(bowtie2_command, stderr=handle_sam_report_genome_global)

    handle_sam_report_genome_global.close()

    # convert sam to bam
    sam_to_bam_genome_global_command = ['samtools', 'view', '-Sb', file_sam_genome_global]

    handle_file_bam_genome_global = open(file_bam_genome_global, 'wb')

    subprocess.call(sam_to_bam_genome_global_command, stdout=handle_file_bam_genome_global)

    # sort bam files
    sort_bam_genome_global_command = ['samtools', 'sort', file_bam_genome_global, '-o', file_sorted_bam_genome_global]

    subprocess.call(sort_bam_genome_global_command)

    # Clean up
    if keep_sam == 0:
        os.remove(file_sam_genome_global)
        print('sam file deleted')

    
    os.remove(file_bam_genome_global)

    # Create bam index files
    create_bam_genome_global_index_command = ['samtools', 'index', file_sorted_bam_genome_global]
    subprocess.call(create_bam_genome_global_index_command)

    os.chdir(initial_dir)



#############################################################

def guideseq_filtered_mapq_AS_primary(dir_sample, global_or_local, amplicon_info, single_or_paired = 'single', min_MAPQ = 25, min_AS=-180):
    #function beginning
    # define samples
    N7 = amplicon_info['index_I1']
    N5 = amplicon_info['index_I2']
    name = amplicon_info['name']


    
    if global_or_local == 'global':
        if single_or_paired == 'single':
            bam_file = create_filename(dir_sample, N7, N5, 'sorted_bam_genome_global_single')
            final_trimmed_bam_filtered_file = create_filename(dir_sample, N7, N5, 'filtered_and_sorted_genome_global_single')
        elif single_or_paired == 'paired':
            bam_file = create_filename(dir_sample, N7, N5, 'sorted_bam_genome_global_paired')
            final_trimmed_bam_filtered_file = create_filename(dir_sample, N7, N5, 'filtered_and_sorted_genome_global_paired')
 
    if global_or_local == 'local':
        if single_or_paired == 'single':
            bam_file = create_filename(dir_sample, N7, N5, 'sorted_bam_genome_local_single')
            final_trimmed_bam_filtered_file = create_filename(dir_sample, N7, N5, 'filtered_and_sorted_genome_local_single')
        elif single_or_paired == 'paired':
            print('wooo')
            bam_file = create_filename(dir_sample, N7, N5, 'sorted_bam_genome_local_paired')
            final_trimmed_bam_filtered_file = create_filename(dir_sample, N7, N5, 'filtered_and_sorted_genome_local_paired')    
    
    all_reads_unfiltered = []
    all_reads_pass_filter = []
    
    bam_alignment_file = pysam.AlignmentFile(bam_file, 'rb')
    bam_in = bam_alignment_file.fetch()
    
    #make output file
    filtered_reads = pysam.AlignmentFile(final_trimmed_bam_filtered_file, "wb", template=bam_alignment_file)
    
    #test each read
    for read in bam_in:
        all_reads_unfiltered.append(read.query_name)
        if read.has_tag('AS'):
            read_AS = read.get_tag('AS')
        #check read quality is good
        if read.mapping_quality >= min_MAPQ and read_AS >= min_AS and not read.is_secondary:
            all_reads_pass_filter.append(read.query_name)
            filtered_reads.write(read)
    filtered_reads.close()
    bam_alignment_file.close()


    count_unfiltered = len(all_reads_unfiltered)

    count_filtered = len(all_reads_pass_filter)
    
    
    
    results_df = pd.DataFrame({'sample_name': [name],
                               'i7': [N7],
                               'i5': [N5],
                               'single_or_paired_reads' : [single_or_paired],
                               'all_alignments_count': [count_unfiltered],
                               'filtered_alignments_count': [count_filtered]
                               },
                                  columns=['sample_name',
                                           'i7',
                                           'i5',
                                           'single_or_paired_reads',
                                           'all_alignments_count',
                                           'filtered_alignments_count'])
    print('single por paired', single_or_paired)
    
    print(results_df)   
    return results_df



###################### TRIMMING for AAV-seq ################## ####
#need to trim off the 3' end of the short reads. This is for amplicons that were too short and have the other side on them.
#      expect more trimmed read 2 and than read 1. Read 2 is longer htan read one as we trim off the ~30 bases of good priming 
#     form read1 so sort amplicons get trimmed for read2 before being trimmed for read1.
#
#    I had wanted to trim 5' to the litated interface, butt this is not easy because the breaksite is variable depending how AAV intgrated    

def trim_AAVseq(dir_sample, amplicon_info, direction5 = 'ATACCGTTATTAACATATGACAACTCAATTAAAC'):

#### need to make this more like my original trim to break sequence end-to-end alignment of replace targeting
    
    # Read1  This is the gene specific primering binding
    R1primer = 'TCTCTGCGCGCTCGCTCGCTCACTGA'  #for the i5 side end of read2 for tn5
    R1_primer_rev = reverse_complement(R1primer)
    
    # reverse trimming This is the Tn5 primer sequence which is acutally the illumina adapter sequence
   # nexteraR1 ='GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG'
    nexteraR1 ='GGGCTCGGAGATGTGTATAAGAGACAG'

    nexteraR1_rev = reverse_complement(nexteraR1)
    
    # We first check if the experiment had any guides
    N7 = amplicon_info['index_I1']
    N5 = amplicon_info['index_I2']
    
    file_R1 = create_filename(dir_sample, N7, N5, 'R1fastq_CorrPrime')
    file_R2 = create_filename(dir_sample, N7, N5, 'R2fastq_CorrPrime')

    file_cutadapt_R1 = create_filename(dir_sample, N7, N5, 'R1trimmed')
    file_cutadapt_R2 = create_filename(dir_sample, N7, N5, 'R2trimmed')
    file_cutadapt_report = create_filename(dir_sample, N7, N5, 'trimmed_report')
    
    
    if not os.path.exists(os.path.dirname(file_cutadapt_R1)):
        os.mkdir(os.path.dirname(file_cutadapt_R1))

    # remove adapters with cutadapt
    #original uditas peramiter had an error -e 0.33 (but was cutting of random stuff too much)
    # -a is hte 3' adapter for Read1
    # -A is 3' adapter for Read2
    # -m minium length
    cutadapt_command = ['cutadapt',
                        '-m', '50',
                        '-e', '0.20',
                        '-a', nexteraR1_rev,
                        '-A', R1_primer_rev,
                        '-o', file_cutadapt_R1, '-p', file_cutadapt_R2,
                        file_R1, file_R2]

    handle_cutadapt_report = open(file_cutadapt_report, 'wb')
    subprocess.call(cutadapt_command, stdout=handle_cutadapt_report)
    handle_cutadapt_report.close()



