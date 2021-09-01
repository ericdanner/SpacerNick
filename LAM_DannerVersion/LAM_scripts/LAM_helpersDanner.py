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
from fastDamerauLevenshtein import damerauLevenshtein


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
        
########################################################################
# Function to make the 'amplicon_info' list. Taking the line of our experiment csv file

def get_csv_data(dir_sample, line_of_data_from_csv):
    sample_info_filename = os.path.join(dir_sample, 'sample_info.csv')
    experiments = pd.read_csv(sample_info_filename)
    return experiments.loc[line_of_data_from_csv]

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
# Select closest barcode with a maximum number of mismatches
# By default it returns barcodes with a maximum of n_max_mismatches=2 mismatches
################################################################################
def select_barcode(seq, barcode_list, n_max_mismatches=1):
    # This compares with all barcodes and selects the one with the smallest hamming distance
    # Before calling this function check if the sequence is already a barcode
    matched_barcodes = list()
    distances = list()
    for barcode in barcode_list:
        h_d = hamm_dist(seq, barcode)
        if h_d <= n_max_mismatches:
            matched_barcodes.append(barcode)
            distances.append(h_d)
    indices = [i for i, x in enumerate(distances) if x == min(distances)]
    return [matched_barcodes[i] for i in indices]


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
# Create umi dict
################################################################################
def create_umi_dict(filename):

    umi_file = open_fastq_or_gz(filename)

    umi_dict = dict()

    umi_reads = itertools.izip(umi_file)

    for header_umi in umi_reads:

        seq_umi = umi_reads.next()
        umi_reads.next()
        qual_umi = umi_reads.next()
        umi_dict[header_umi[0].split()[0][1:]] = [seq_umi[0].rstrip(), qual_umi[0].rstrip()]

    return umi_dict



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
    elif filetype == 'R1fastq':
        return os.path.join(main_folder, 'fastq_files', N7 + '_' + N5 + '_R1.fastq')
    elif filetype == 'R1fastqgz':
        return os.path.join(main_folder, 'fastq_files', N7 + '_' + N5 + '_R1.fastq.gz')
    elif filetype == 'R2fastq':
        return os.path.join(main_folder, 'fastq_files', N7 + '_' + N5 + '_R2.fastq')
    elif filetype == 'R2fastqgz':
        return os.path.join(main_folder, 'fastq_files', N7 + '_' + N5 + '_R2.fastq.gz')
   
    ### Added to make LAM cleanup after extracting UMI and cleaving adapter
    elif filetype == 'R1fastq_LAM':
        return os.path.join(main_folder, 'fastq_files', N7 + '_' + N5 + '_R1.LAM.fastq')
    elif filetype == 'R1fastqgz_LAM':
        return os.path.join(main_folder, 'fastq_files', N7 + '_' + N5 + '_R1.LAM.fastq.gz')
    elif filetype == 'R2fastq_LAM':
        return os.path.join(main_folder, 'fastq_files', N7 + '_' + N5 + '_R2.LAM.fastq')
    elif filetype == 'R2fastqgz_LAM':
        return os.path.join(main_folder, 'fastq_files', N7 + '_' + N5 + '_R2.LAM.fastq.gz')
    
    elif filetype == 'breaktrimmed_R1fastq':
        return os.path.join(main_folder, 'fastq_files', N7 + '_' + N5 + '_R1.breaktrimmed.fastq')
    elif filetype == 'breaktrimmed_R1fastqgz':
        return os.path.join(main_folder, 'fastq_files', N7 + '_' + N5 + '_R1.breaktrimmed.fastq.gz')
    elif filetype == 'breaktrimmed_R2fastq':
        return os.path.join(main_folder, 'fastq_files', N7 + '_' + N5 + '_R2.breaktrimmed.fastq')
    elif filetype == 'breaktrimmed_R2fastqgz':
        return os.path.join(main_folder, 'fastq_files', N7 + '_' + N5 + '_R2.breaktrimmed.fastq.gz')
    
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
    #####
    
    elif filetype == 'umifastq':
        return os.path.join(main_folder, 'fastq_files', N7 + '_' + N5 + '_umi.fastq')
    elif filetype == 'umifastqgz':
        return os.path.join(main_folder, 'fastq_files', N7 + '_' + N5 + '_umi.fastq.gz')
    
    ### Added these for trimming to breaksite
    elif filetype == 'R1trimmed2break':
        return os.path.join(main_folder, 'cutadapt_files_to_breaksite', N7 + '_' + N5 + '_R1.trimmed2break.fastq.gz')
    elif filetype == 'R2trimmed2break':
        return os.path.join(main_folder, 'cutadapt_files_to_breaksite', N7 + '_' + N5 + '_R2.trimmed2break.fastq.gz')
    elif filetype == 'trimmed_report2break':
        return os.path.join(main_folder, 'cutadapt_files_to_breaksite', N7 + '_' + N5 + '.trimmed.report.txt')
    # added. for aligning to amplicons first
    elif filetype == 'break_trimmed_sam_amplicons':         
        return os.path.join(main_folder, 'break_trimmed_sam_amplicon_files', N7 + '_' + N5 + '.sam')
    elif filetype == 'break_trimmed_sam_report_amplicons':
        return os.path.join(main_folder, 'break_trimmed_sam_amplicon_files', N7 + '_' + N5 + '.sam.report.txt')
    elif filetype == 'break_trimmed_bam_amplicons':
        return os.path.join(main_folder, 'break_trimmed_bam_amplicon_files', N7 + '_' + N5 + '.bam')
    elif filetype == 'break_trimmed_sorted_bam_amplicons':
        return os.path.join(main_folder, 'break_trimmed_bam_amplicon_files', N7 + '_' + N5 + '.sorted.bam')
    elif filetype == 'break_trimmed_sorted_bai_amplicons':
        return os.path.join(main_folder, 'break_trimmed_bam_amplicon_files', N7 + '_' + N5 + '.sorted.bam.bai')
    elif filetype == 'break_trimmed_unmapped_bam_amplicons':
        return os.path.join(main_folder, 'break_trimmed_bam_amplicon_files', N7 + '_' + N5 + '_amplicons_unmapped.bam')
    elif filetype == 'break_trimmed_qsorted_unmapped_bam_amplicons':
        return os.path.join(main_folder, 'break_trimmed_bam_amplicon_files', N7 + '_' + N5 + '_qsorted_amplicons_unmapped.bam')
    elif filetype == 'break_trimmed_unmapped_amplicons_R1fastq':
        return os.path.join(main_folder, 'break_trimmed_amplicons_unmapped_fastq_files', N7 + '_' + N5 + '_amplicons_unmapped_R1.fastq')
    elif filetype == 'break_trimmed_unmapped_amplicons_R2fastq':
        return os.path.join(main_folder, 'break_trimmed_amplicons_unmapped_fastq_files', N7 + '_' + N5 + '_amplicons_unmapped_R2.fastq')
    elif filetype == 'break_trimmed_unmapped_amplicons_R1fastqgz':
        return os.path.join(main_folder, 'break_trimmed_amplicons_unmapped_fastq_files', N7 + '_' + N5 + '_amplicons_unmapped_R1.fastq.gz')
    elif filetype == 'break_trimmed_unmapped_amplicons_R2fastqgz':
        return os.path.join(main_folder, 'break_trimmed_amplicons_unmapped_fastq_files', N7 + '_' + N5 + '_amplicons_unmapped_R2.fastq.gz')
    elif filetype == 'break_trimmed_unmapped_amplicons_report':
        return os.path.join(main_folder, 'break_trimmed_amplicons_unmapped_fastq_files', N7 + '_' + N5 + '.unmapped.report.txt')
    elif filetype == 'break_trimmed_sam_genome_global':
        return os.path.join(main_folder, 'break_trimmed_sam_genome_global_files', N7 + '_' + N5 + '.sam')
    elif filetype == 'break_trimmed_sam_report_genome_global':
        return os.path.join(main_folder, 'break_trimmed_sam_genome_global_files', N7 + '_' + N5 + '.sam.report.txt')
    elif filetype == 'break_trimmed_bam_genome_global':
        return os.path.join(main_folder, 'break_trimmed_bam_genome_global_files', N7 + '_' + N5 + '.bam')
    elif filetype == 'break_trimmed_sorted_bam_genome_global':
        return os.path.join(main_folder, 'break_trimmed_bam_genome_global_files', N7 + '_' + N5 + '.sorted.bam')
    elif filetype == 'break_trimmed_sorted_bai_genome_global':
        return os.path.join(main_folder, 'break_trimmed_bam_genome_global_files', N7 + '_' + N5 + '.sorted.bam.bai')
    elif filetype == 'break_trimmed_filtered_and_sorted_genome_global':
        return os.path.join(main_folder, 'break_trimmed_bam_genome_global_files', N7 + '_' + N5 + '_primary_as_mapq.sorted.bam')    
    elif filetype == 'break_trimmed_genome_global_bed':
        return os.path.join(main_folder, 'break_trimmed_bam_genome_global_files', N7 + '_' + N5 + 'global.sorted.bed')
    elif filetype == 'break_trimmed_sam_genome_local':
        return os.path.join(main_folder, 'break_trimmed_sam_genome_local_files', N7 + '_' + N5 + '.sam')
    elif filetype == 'break_trimmed_sam_report_genome_local':
        return os.path.join(main_folder, 'break_trimmed_sam_genome_local_files', N7 + '_' + N5 + '.sam.report.txt')
    elif filetype == 'break_trimmed_bam_genome_local':
        return os.path.join(main_folder, 'break_trimmed_bam_genome_local_files', N7 + '_' + N5 + '.bam')
    elif filetype == 'break_trimmed_sorted_bam_genome_local':
        return os.path.join(main_folder, 'break_trimmed_bam_genome_local_files', N7 + '_' + N5 + '.sorted.bam')
    elif filetype == 'break_trimmed_sorted_bai_genome_local':
        return os.path.join(main_folder, 'break_trimmed_bam_genome_local_files', N7 + '_' + N5 + '.sorted.bam.bai')
    elif filetype == 'break_trimmed_filtered_and_sorted_bam_genome_local':
        return os.path.join(main_folder, 'break_trimmed_bam_genome_local_files', N7 + '_' + N5 + '_primary_as_mapq.sorted.bam')
    elif filetype == 'break_trimmed_genome_local_bed':
        return os.path.join(main_folder, 'break_trimmed_bam_genome_local_files', N7 + '_' + N5 + 'local.sorted.bed')
    elif filetype == 'all_bed':
        return os.path.join(main_folder, '..', 'bed_files')
    ###############################################
    
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
    elif filetype == 'sam_AAV_local':
        return os.path.join(main_folder, 'sam_AAV_local_files', N7 + '_' + N5 + '.sam')
    elif filetype == 'sam_report_AAV_local':
        return os.path.join(main_folder, 'sam_AAV_local_files', N7 + '_' + N5 + '.sam.report.txt')
    elif filetype == 'bam_plasmid_local':
        return os.path.join(main_folder, 'bam_plasmid_local_files', N7 + '_' + N5 + '.bam')
    elif filetype == 'sorted_bam_plasmid_local':
        return os.path.join(main_folder, 'bam_plasmid_local_files', N7 + '_' + N5 + '.sorted.bam')
    elif filetype == 'bam_AAV_local':
        return os.path.join(main_folder, 'bam_AAV_local_files', N7 + '_' + N5 + '.bam')
    elif filetype == 'sorted_bam_AAV_local':
        return os.path.join(main_folder, 'bam_AAV_local_files', N7 + '_' + N5 + '.sorted.bam')
    elif filetype == 'sorted_bai_plasmid_local':
        return os.path.join(main_folder, 'bam_plasmid_local_files', N7 + '_' + N5 + '.sorted.bam.bai')
    elif filetype == 'unmapped_bam_plasmid_local':
        return os.path.join(main_folder, 'bam_plasmid_local_files', N7 + '_' + N5 + '_unmapped.bam')
    elif filetype == 'mapped_bam_plasmid_local':
        return os.path.join(main_folder, 'bam_plasmid_local_files', N7 + '_' + N5 + '_mapped.bam')
    elif filetype == 'qsorted_unmapped_bam_plasmid_local':
        return os.path.join(main_folder, 'bam_plasmid_local_files', N7 + '_' + N5 + '_qsorted_unmapped.bam')
    elif filetype == 'qsorted_mapped_bam_plasmid_local':
        return os.path.join(main_folder, 'bam_plasmid_local_files', N7 + '_' + N5 + '_qsorted_mapped.bam')
    elif filetype == 'unmapped_plasmid_R1fastq':
        return os.path.join(main_folder, 'plasmid_unmapped_fastq_files', N7 + '_' + N5 + '_plasmid_unmapped_R1.fastq')
    elif filetype == 'unmapped_plasmid_R2fastq':
        return os.path.join(main_folder, 'plasmid_unmapped_fastq_files', N7 + '_' + N5 + '_plasmid_unmapped_R2.fastq')
    elif filetype == 'mapped_plasmid_R1fastq':
        return os.path.join(main_folder, 'plasmid_mapped_fastq_files', N7 + '_' + N5 + '_plasmid_mapped_R1.fastq')
    elif filetype == 'mapped_plasmid_R2fastq':
        return os.path.join(main_folder, 'plasmid_mapped_fastq_files', N7 + '_' + N5 + '_plasmid_mapped_R2.fastq')
    elif filetype == 'unmapped_plasmid_R1fastqgz':
        return os.path.join(main_folder, 'plasmid_unmapped_fastq_files', N7 + '_' + N5 + '_plasmid_unmapped_R1.fastq.gz')
    elif filetype == 'unmapped_plasmid_R2fastqgz':
        return os.path.join(main_folder, 'plasmid_unmapped_fastq_files', N7 + '_' + N5 + '_plasmid_unmapped_R2.fastq.gz')
    elif filetype == 'mapped_plasmid_R1fastqgz':
        return os.path.join(main_folder, 'plasmid_mapped_fastq_files', N7 + '_' + N5 + '_plasmid_mapped_R1.fastq.gz')
    elif filetype == 'mapped_plasmid_R2fastqgz':
        return os.path.join(main_folder, 'plasmid_mapped_fastq_files', N7 + '_' + N5 + '_plasmid_mapped_R2.fastq.gz')
    
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
    
###########################
#
#   my homemade function to make a fastq with only correct priming events
#   I currently have skipped making this becasue a relitivly high amount are on target due to nested and I just want results
#   this should be a similiar function to the demultiplex sorting function but simpler
#   the sequence to align is the first reads of r2
#######################

def correct_priming(dir_sample, amplicon_info, primer_seq_plus_downstream):
    
    #defined the sequence as input (primer and 12 nt downstream). Normally 32-37nt sequence
    n_max_mismatches = 3
    length_primer_down = len(primer_seq_plus_downstream)
    
    files_out = list()
    n_file = 0
    files_out_dict = dict()

    for file_selected in files_out:
        files_out_dict[os.path.basename(file_selected)] = n_file
        n_file += 1
    
    
    #define the input files and output files
    N7 = amplicon_info['index_I1']
    N5 = amplicon_info['index_I2']
    
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

        reads_in_experiment_list_count = 0

        reads_misprimed = 0

        reads_good_prime = 0
    
 
        r1_r2 = itertools.izip(r1_file, r2_file)

        for header_r1, header_r2 in r1_r2:
            reads_in_experiment_list_count += 1
            
    
            seq_r1, seq_r2 = r1_r2.next()

            r1_r2.next()

            qual_r1, qual_r2 = r1_r2.next()

            seq_r1, seq_r2 = seq_r1.rstrip(), seq_r2.rstrip()

            qual_r1, qual_r2 = qual_r1.rstrip(), qual_r2.rstrip()

            #We mask with N any bases with scores below or equal to , (11, default in mask)
            
            #don't need seq1
            #seq_r1_use = mask(seq_r1, qual_r1)

            #this is modified for miniseq UMI->Index2
            seq_r2_use = mask(seq_r2[:length_primer_down], qual_r2[:length_primer_down])

    
            # change to 1 for reads with perfect match or match within hamming distance decided above
            is_good_read = 0

            if seq_r2_use == primer_seq_plus_downstream:
                # perfect match case
                is_good_read = 1
            else:
                # We look for barcodes with up to two mismatches, default in select_barcode
                h_d = hamm_dist(seq_r2_use, primer_seq_plus_downstream)
                if h_d <= n_max_mismatches:
                    # match after checking within hamming distance
                    is_good_read = 1
            
                        
            #Print good reads 
            if is_good_read:
                reads_good_prime += 1
                
                # We test whether the read has on of the combination of indices from our experiment list
                # If not save in a separate file
                r1f = ref_file_R1_corrprime
                r2f = ref_file_R2_corrprime

                print("\n".join([header_r1.rstrip(), seq_r1.rstrip(), "+", qual_r1.rstrip()]), file=r1f)
                print("\n".join([header_r2.rstrip(), seq_r2.rstrip(), "+", qual_r2.rstrip()]), file=r2f)
                    
                
                

            else:
                # We print reads with mispriming
                print("\n".join([header_r1.rstrip(), seq_r1.rstrip(), "+", qual_r1.rstrip()]), file=ref_file_out_r1_not_corrprime)
                print("\n".join([header_r2.rstrip(), seq_r2.rstrip(), "+", qual_r2.rstrip()]), file=ref_file_out_r2_not_corrprime)
                
                reads_misprimed += 1
                
        print('reads_in_experiment_list_count', reads_in_experiment_list_count)
        print(' reads with good priming:', reads_good_prime)
        print('reads misprimed', reads_misprimed)


   ############################
    #
    #  Function to determine the kind of reaction from the number of cuts and their locations
    # I changed this to add the replace option
    #
    # The function classify the cases:
    #   - No cut (just UDiTaS primer, used for controls)
    #   - Single cut
    #   - Replace (has to come before duel cut in logic order)
    #   - Dual cut on same chromosome, generates amplicons with large deletions, etc
    #      > Duel cut on same chromosome can also generate HDR product.
    #   - Dual cut on different chromosomes. Generates 10 amplicons including translocations
    #   - Triple cuts on different chromosomes. Generates 21 amplicons including translocations. NOTE that if two of the
    #     cuts are in the same chromosome and close together (less than amplicon_window_around_cut) the results may
    #     be incorrect since some reads may be mapped to multiple amplicons, but only counted around the cut in one amplicon
    #
    ############################
def get_reaction_type(amplicon_info):
    has_guide1 = type(amplicon_info['chr_guide_1']) is str or type(amplicon_info['chr_guide_1']) is unicode
    has_guide2 = type(amplicon_info['chr_guide_2']) is str or type(amplicon_info['chr_guide_2']) is unicode
    has_guide3 = type(amplicon_info['chr_guide_3']) is str or type(amplicon_info['chr_guide_3']) is unicode
    has_replace_donor = type(amplicon_info['Replace_Donor']) is str or type(amplicon_info['Replace_Donor']) is unicode
    has_HDR = type(amplicon_info['HDR_correct_seq']) is str or type(amplicon_info['HDR_correct_seq']) is unicode

    
    if not has_guide1 and not has_guide2 and not has_guide3:
        reaction_type = 'control'
    elif has_guide1 and not has_guide2 and not has_guide3:
        reaction_type = 'single_cut'
     #this is the modification of the origianl. It shows if we are doing a replace targeting
    elif has_guide1 and has_guide2 and amplicon_info['chr_guide_1'] == amplicon_info['chr_guide_2'] and has_replace_donor:
        reaction_type = 'replace' 
    
    elif has_guide1 and has_guide2 and amplicon_info['chr_guide_1'] == amplicon_info['chr_guide_2'] and not has_guide3:
        if has_HDR:
            reaction_type = 'double_cut_same_chromosome_and_HDR'
        else:
            reaction_type = 'double_cut_same_chromosome'

    elif has_guide1 and has_guide2 and amplicon_info['chr_guide_1'] != amplicon_info['chr_guide_2'] and not has_guide3:
        reaction_type = 'double_cut_different_chromosomes'
    elif has_guide1 and has_guide2 and has_guide3:
        if (amplicon_info['chr_guide_1'] == amplicon_info['chr_guide_2'] or
                amplicon_info['chr_guide_1'] == amplicon_info['chr_guide_3'] or
                amplicon_info['chr_guide_2'] == amplicon_info['chr_guide_3']):
            raise ReactionTypeError('The reaction with three cuts with at least two in the same chromosome is' +
                                    ' not yet supported by current version of UDiTaS')
        reaction_type = 'triple_cut_different_chromosomes'
    else:
        raise ReactionTypeError('Reaction type not yet supported by current version of UDiTaS')

    return reaction_type
    
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

################################################
# Function to write reference plasmid sequence
################################################
def create_plasmid_reference(dir_sample, amplicon_info):
    N7 = amplicon_info['index_I1']
    N5 = amplicon_info['index_I2']
    exp_dir = create_filename(dir_sample, N7, N5, 'mainfolder')
    amplicon_folder = os.path.join(exp_dir, 'amplicons')
    if not os.path.exists(amplicon_folder):
        os.mkdir(amplicon_folder)

    filename = os.path.join(exp_dir, amplicon_folder, 'plasmid.fa')
    file_handle = open(filename, "w")

    # If several plasmids were given, separate them using ';' in sample_info.csv. Here we remove the ';' so that
    # we just concatenate the sequences
    pl_seq = amplicon_info['plasmid_sequence'].replace(';', '')
    seq1 = Seq(pl_seq, IUPAC.unambiguous_dna)
    record1 = SeqRecord(seq1, 'plasmid', description='')
    SeqIO.write(record1, file_handle, 'fasta')

    file_handle.close()
    # Create index file
    initial_dir = os.getcwd()
    os.chdir(amplicon_folder)
    index_err_file = os.path.join(amplicon_folder, 'index_plasmid.err')
    index_out_file = os.path.join(amplicon_folder, 'index_plasmid.out')

    index_err_fh = open(index_err_file, 'wb')
    index_out_fh = open(index_out_file, 'wb')
    subprocess.call(['bowtie2-build',
                     filename, 'plasmid'], stderr=index_err_fh, stdout=index_out_fh)
    os.chdir(initial_dir)
    index_err_fh.close()
    index_out_fh.close() 
    

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


############################
#
# Create amplicon. Creates fasta file with the custom reference amplicons including deletions, inversions, etc...
#
# Input: dir_sample, directory to be analyzed
#        amplicon_info, slice of sample_info.csv for the sample being processed
#        file_genome_2bit, 2bit file with the reference genome being used
#        amplicon_window_around_cut, value used to grab sequences around cut sites
#
#   This functions creates a set of reference amplicons built from the expected fragments after cutting
#   The amplicons are built differently depending on the case classified by get_reaction_type()
#
# ##########################
def create_amplicon(dir_sample, amplicon_info, file_genome_2bit, amplicon_window_around_cut=1000):
    # We first check the reaction type
    reaction_type = get_reaction_type(amplicon_info)

    genome = twobitreader.TwoBitFile(file_genome_2bit)  # Load genome. Used for getting the sequences

    amplicon_list = []

    # For all reaction types, we check that we don't go out of boundaries when building the amplicons
    # This is unlikely for hg38 or mm10, but could easily happen in the UDiTaS primer is in a plasmid
    if reaction_type == 'control':
        # Case no guides
        if amplicon_info['strand'] == '+':  # This is the UDiTaS oligo strand
            end_coordinate = int(amplicon_info['start']) + amplicon_window_around_cut
            if end_coordinate > len(genome[amplicon_info['chr']]):
                end_coordinate = len(genome[amplicon_info['chr']])
            amplicon_list.append(['wt', genome[amplicon_info['chr']][int(amplicon_info['start']):end_coordinate]])
        elif amplicon_info['strand'] == '-':
            start_coordinate = int(amplicon_info['end']) - amplicon_window_around_cut
            if start_coordinate < 0:
                start_coordinate = 0
            amplicon_list.append(['wt', genome[amplicon_info['chr']][start_coordinate:int(amplicon_info['end'])]])
        else:
            raise StrandError('strand can only have as values + or -')
    elif reaction_type == 'single_cut':
        # Case one guide
        if amplicon_info['strand_guide_1'] == '+':
            # sp or sa for the moment only
            cut1 = amplicon_info['end_guide_1'] - 3
        elif amplicon_info['strand_guide_1'] == '-':
            cut1 = amplicon_info['start_guide_1'] + 3
        else:
            raise StrandError('strand_guide_1 can only have as values + or -')

        start_coordinate = int(cut1 - amplicon_window_around_cut)
        if start_coordinate < 0:
            start_coordinate = 0
        end_coordinate = int(cut1 + amplicon_window_around_cut)
        if end_coordinate > len(genome[amplicon_info['chr_guide_1']]):
            end_coordinate = len(genome[amplicon_info['chr_guide_1']])

        seq_upstream = genome[amplicon_info['chr_guide_1']][start_coordinate:int(cut1)]
        seq_downstream = genome[amplicon_info['chr_guide_1']][int(cut1):end_coordinate]

        amplicon_list.append(['wt', seq_upstream + seq_downstream])
        amplicon_list.append(['1a_1a', seq_upstream + reverse_complement(seq_upstream)])
        amplicon_list.append(['1b_1b', reverse_complement(seq_downstream) + seq_downstream])
    elif reaction_type == 'double_cut_same_chromosome':
        # Case two guides on the same chromosome
        if amplicon_info['strand_guide_1'] == '+':
            # sp or sa for the moment only
            cut1 = amplicon_info['end_guide_1'] - 3
        elif amplicon_info['strand_guide_1'] == '-':
            cut1 = amplicon_info['start_guide_1'] + 3
        else:
            raise StrandError('strand_guide_1 can only have as values + or -')

        if amplicon_info['strand_guide_2'] == '+':
            cut2 = amplicon_info['end_guide_2'] - 3
        elif amplicon_info['strand_guide_2'] == '-':
            cut2 = amplicon_info['start_guide_2'] + 3
        else:
            raise StrandError('strand_guide_2 can only have as values + or -')

        # We switch the coordinates of cut1 and cut2 if the guides are provided so that cut2 < cut1
        # cut1 it will always be smaller than cut2 and in the results cut1 will be the cut site with
        # smaller genomic coordinate
        # cut1 and cut2 also flipped in get_cut_in_reference_amplicon_df

        if cut2 < cut1:
            (cut1, cut2) = (cut2, cut1)

        start_coordinate = int(cut1 - amplicon_window_around_cut)
        if start_coordinate < 0:
            start_coordinate = 0
        end_coordinate = int(cut2 + amplicon_window_around_cut)
        if end_coordinate > len(genome[amplicon_info['chr_guide_1']]):
            end_coordinate = len(genome[amplicon_info['chr_guide_1']])

        seq_upstream = genome[amplicon_info['chr_guide_1']][start_coordinate:int(cut1)]
        seq_cut1_cut2 = genome[amplicon_info['chr_guide_1']][int(cut1):int(cut2)]
        seq_downstream = genome[amplicon_info['chr_guide_1']][int(cut2):end_coordinate]

        amplicon_list.append(['wt', seq_upstream + seq_cut1_cut2 + seq_downstream])
        amplicon_list.append(['large_deletion', seq_upstream + seq_downstream])
        amplicon_list.append(['large_inversion', seq_upstream + reverse_complement(seq_cut1_cut2) + seq_downstream])
        amplicon_list.append(['1a_1a', seq_upstream + reverse_complement(seq_upstream)])
        amplicon_list.append(['2b_2b', reverse_complement(seq_downstream) + seq_downstream])
        
    elif reaction_type == 'double_cut_same_chromosome_and_HDR':
        
        # prepare the HDR sequence
        HDR_seq = amplicon_info['HDR_correct_seq'].replace(';', '')
        
        # Case two guides on the same chromosome
        if amplicon_info['strand_guide_1'] == '+':
            # sp or sa for the moment only
            cut1 = amplicon_info['end_guide_1'] - 3
        elif amplicon_info['strand_guide_1'] == '-':
            cut1 = amplicon_info['start_guide_1'] + 3
        else:
            raise StrandError('strand_guide_1 can only have as values + or -')

        if amplicon_info['strand_guide_2'] == '+':
            cut2 = amplicon_info['end_guide_2'] - 3
        elif amplicon_info['strand_guide_2'] == '-':
            cut2 = amplicon_info['start_guide_2'] + 3
        else:
            raise StrandError('strand_guide_2 can only have as values + or -')

        # We switch the coordinates of cut1 and cut2 if the guides are provided so that cut2 < cut1
        # cut1 it will always be smaller than cut2 and in the results cut1 will be the cut site with
        # smaller genomic coordinate
        # cut1 and cut2 also flipped in get_cut_in_reference_amplicon_df

        if cut2 < cut1:
            (cut1, cut2) = (cut2, cut1)

        start_coordinate = int(cut1 - amplicon_window_around_cut)
        if start_coordinate < 0:
            start_coordinate = 0
        end_coordinate = int(cut2 + amplicon_window_around_cut)
        if end_coordinate > len(genome[amplicon_info['chr_guide_1']]):
            end_coordinate = len(genome[amplicon_info['chr_guide_1']])

        seq_upstream = genome[amplicon_info['chr_guide_1']][start_coordinate:int(cut1)]
        seq_cut1_cut2 = genome[amplicon_info['chr_guide_1']][int(cut1):int(cut2)]
        seq_downstream = genome[amplicon_info['chr_guide_1']][int(cut2):end_coordinate]

        amplicon_list.append(['wt', seq_upstream + seq_cut1_cut2 + seq_downstream])
        amplicon_list.append(['large_deletion', seq_upstream + seq_downstream])
        amplicon_list.append(['large_inversion', seq_upstream + reverse_complement(seq_cut1_cut2) + seq_downstream])
        amplicon_list.append(['1a_1a', seq_upstream + reverse_complement(seq_upstream)])
        amplicon_list.append(['2b_2b', reverse_complement(seq_downstream) + seq_downstream])
        amplicon_list.append(['HDR_knockin', HDR_seq])
    
    #added this modification from Uditas for replace
    elif reaction_type == 'replace':
        if amplicon_info['strand_guide_1'] == '+':
            # sp or sa for the moment only2b
            cut1 = amplicon_info['end_guide_1'] - 3
        elif amplicon_info['strand_guide_1'] == '-':
            cut1 = amplicon_info['start_guide_1'] + 3
        else:
            raise StrandError('strand_guide_1 can only have as values + or -')

        if amplicon_info['strand_guide_2'] == '+':
            cut2 = amplicon_info['end_guide_2'] - 3
        elif amplicon_info['strand_guide_2'] == '-':
            cut2 = amplicon_info['start_guide_2'] + 3
        else:
            raise StrandError('strand_guide_2 can only have as values + or -')

        # We switch the coordinates of cut1 and cut2 if the guides are provided so that cut2 < cut1
        # cut1 it will always be smaller than cut2 and in the results cut1 will be the cut site with
        # smaller genomic coordinate
        # cut1 and cut2 also flipped in get_cut_in_reference_amplicon_df

        if cut2 < cut1:
            (cut1, cut2) = (cut2, cut1)

        start_coordinate = int(cut1 - amplicon_window_around_cut)
        if start_coordinate < 0:
            start_coordinate = 0
        end_coordinate = int(cut2 + amplicon_window_around_cut)
        if end_coordinate > len(genome[amplicon_info['chr_guide_1']]):
            end_coordinate = len(genome[amplicon_info['chr_guide_1']])

        seq_upstream = genome[amplicon_info['chr_guide_1']][start_coordinate:int(cut1)]
        seq_cut1_cut2 = genome[amplicon_info['chr_guide_1']][int(cut1):int(cut2)]
        seq_downstream = genome[amplicon_info['chr_guide_1']][int(cut2):end_coordinate]

#this does not cover everything because for hiti if you invert the excised exon it gets confusing so I jus tlook at the replace-seq + excised exon junction
        amplicon_list.append(['wt', seq_upstream + seq_cut1_cut2 + seq_downstream])
        amplicon_list.append(['large_deletion', seq_upstream + seq_downstream])
        amplicon_list.append(['large_inversion', seq_upstream + reverse_complement(seq_cut1_cut2) + seq_downstream])
        amplicon_list.append(['replace_fwd', seq_upstream + amplicon_info['Replace_Donor'] + seq_downstream])
        amplicon_list.append(['replace_rev', seq_upstream + reverse_complement(amplicon_info['Replace_Donor']) + seq_downstream])
        amplicon_list.append(['hiti_5prime_replace_exon_fwd',  amplicon_info['Replace_Donor'] + seq_cut1_cut2])
        amplicon_list.append(['hiti_5prime_replace_exon_rev',  reverse_complement(amplicon_info['Replace_Donor']) + seq_cut1_cut2])
        amplicon_list.append(['hiti_3prime_replace_exon_fwd', seq_cut1_cut2 + amplicon_info['Replace_Donor']])
        amplicon_list.append(['hiti_3prime_replace_exon_rev', seq_cut1_cut2 + reverse_complement(amplicon_info['Replace_Donor'])])
        amplicon_list.append(['doner_tail_tail', amplicon_info['Replace_Donor'] + reverse_complement(amplicon_info['Replace_Donor'])])
        amplicon_list.append(['doner_head_tail', amplicon_info['Replace_Donor'] + amplicon_info['Replace_Donor']])
        amplicon_list.append(['doner_head_head', reverse_complement(amplicon_info['Replace_Donor']) + amplicon_info['Replace_Donor']])      
        amplicon_list.append(['1a_1a', seq_upstream + reverse_complement(seq_upstream)])
        amplicon_list.append(['2b_2b', reverse_complement(seq_downstream) + seq_downstream])
        
        
    elif reaction_type == 'double_cut_different_chromosomes':
        # Case two guides on different chromosomes
        if amplicon_info['strand_guide_1'] == '+':
            # sp or sa for the moment only
            cut1 = amplicon_info['end_guide_1'] - 3
        elif amplicon_info['strand_guide_1'] == '-':
            cut1 = amplicon_info['start_guide_1'] + 3
        else:
            raise StrandError('strand_guide_1 can only have as values + or -')

        if amplicon_info['strand_guide_2'] == '+':
            cut2 = amplicon_info['end_guide_2'] - 3
        elif amplicon_info['strand_guide_2'] == '-':
            cut2 = amplicon_info['start_guide_2'] + 3
        else:
            raise StrandError('strand_guide_2 can only have as values + or -')

        start_coordinate1 = int(cut1 - amplicon_window_around_cut)
        if start_coordinate1 < 0:
            start_coordinate1 = 0
        end_coordinate1 = int(cut1 + amplicon_window_around_cut)
        if end_coordinate1 > len(genome[amplicon_info['chr_guide_1']]):
            end_coordinate1 = len(genome[amplicon_info['chr_guide_1']])

        start_coordinate2 = int(cut2 - amplicon_window_around_cut)
        if start_coordinate2 < 0:
            start_coordinate2 = 0
        end_coordinate2 = int(cut2 + amplicon_window_around_cut)
        if end_coordinate2 > len(genome[amplicon_info['chr_guide_2']]):
            end_coordinate2 = len(genome[amplicon_info['chr_guide_2']])

        seq_1a = genome[amplicon_info['chr_guide_1']][start_coordinate1:int(cut1)]
        seq_1b = genome[amplicon_info['chr_guide_1']][int(cut1):end_coordinate1]
        seq_2a = genome[amplicon_info['chr_guide_2']][start_coordinate2:int(cut2)]
        seq_2b = genome[amplicon_info['chr_guide_2']][int(cut2):end_coordinate2]

        amplicon_list.append(['1a_1a', seq_1a + reverse_complement(seq_1a)])
        amplicon_list.append(['1a_1b', seq_1a + seq_1b])
        amplicon_list.append(['1a_2a', seq_1a + reverse_complement(seq_2a)])
        amplicon_list.append(['1a_2b', seq_1a + seq_2b])

        amplicon_list.append(['1b_1b', reverse_complement(seq_1b) + seq_1b])
        amplicon_list.append(['2a_1b', seq_2a + seq_1b])
        amplicon_list.append(['2b_1b', reverse_complement(seq_2b) + seq_1b])

        amplicon_list.append(['2a_2a', seq_2a + reverse_complement(seq_2a)])
        amplicon_list.append(['2a_2b', seq_2a + seq_2b])

        amplicon_list.append(['2b_2b', reverse_complement(seq_2b) + seq_2b])
    elif reaction_type == 'triple_cut_different_chromosomes':
        # Case three guides on different chromosomes

        if amplicon_info['strand_guide_1'] == '+':
            # sp or sa for the moment only
            cut1 = amplicon_info['end_guide_1'] - 3
        elif amplicon_info['strand_guide_1'] == '-':
            cut1 = amplicon_info['start_guide_1'] + 3
        else:
            raise StrandError('strand_guide_1 can only have as values + or -')

        if amplicon_info['strand_guide_2'] == '+':
            cut2 = amplicon_info['end_guide_2'] - 3
        elif amplicon_info['strand_guide_2'] == '-':
            cut2 = amplicon_info['start_guide_2'] + 3
        else:
            raise StrandError('strand_guide_2 can only have as values + or -')

        if amplicon_info['strand_guide_3'] == '+':
            cut3 = amplicon_info['end_guide_3'] - 3
        elif amplicon_info['strand_guide_3'] == '-':
            cut3 = amplicon_info['start_guide_3'] + 3
        else:
            raise StrandError('strand_guide_3 can only have as values + or -')

        start_coordinate1 = int(cut1 - amplicon_window_around_cut)
        if start_coordinate1 < 0:
            start_coordinate1 = 0
        end_coordinate1 = int(cut1 + amplicon_window_around_cut)
        if end_coordinate1 > len(genome[amplicon_info['chr_guide_1']]):
            end_coordinate1 = len(genome[amplicon_info['chr_guide_1']])

        start_coordinate2 = int(cut2 - amplicon_window_around_cut)
        if start_coordinate2 < 0:
            start_coordinate2 = 0
        end_coordinate2 = int(cut2 + amplicon_window_around_cut)
        if end_coordinate2 > len(genome[amplicon_info['chr_guide_2']]):
            end_coordinate2 = len(genome[amplicon_info['chr_guide_2']])

        start_coordinate3 = int(cut3 - amplicon_window_around_cut)
        if start_coordinate3 < 0:
            start_coordinate3 = 0
        end_coordinate3 = int(cut3 + amplicon_window_around_cut)
        if end_coordinate3 > len(genome[amplicon_info['chr_guide_3']]):
            end_coordinate3 = len(genome[amplicon_info['chr_guide_3']])

        seq_1a = genome[amplicon_info['chr_guide_1']][start_coordinate1:int(cut1)]
        seq_1b = genome[amplicon_info['chr_guide_1']][int(cut1):end_coordinate1]
        seq_2a = genome[amplicon_info['chr_guide_2']][start_coordinate2:int(cut2)]
        seq_2b = genome[amplicon_info['chr_guide_2']][int(cut2):end_coordinate2]
        seq_3a = genome[amplicon_info['chr_guide_3']][start_coordinate3:int(cut3)]
        seq_3b = genome[amplicon_info['chr_guide_3']][int(cut3):end_coordinate3]

        amplicon_list.append(['1a_1a', seq_1a + reverse_complement(seq_1a)])
        amplicon_list.append(['1a_1b', seq_1a + seq_1b])
        amplicon_list.append(['1a_2a', seq_1a + reverse_complement(seq_2a)])
        amplicon_list.append(['1a_2b', seq_1a + seq_2b])
        amplicon_list.append(['1a_3a', seq_1a + reverse_complement(seq_3a)])
        amplicon_list.append(['1a_3b', seq_1a + seq_3b])

        amplicon_list.append(['1b_1b', reverse_complement(seq_1b) + seq_1b])
        amplicon_list.append(['2a_1b', seq_2a + seq_1b])
        amplicon_list.append(['2b_1b', reverse_complement(seq_2b) + seq_1b])
        amplicon_list.append(['3a_1b', seq_3a + seq_1b])
        amplicon_list.append(['3b_1b', reverse_complement(seq_3b) + seq_1b])

        amplicon_list.append(['2a_2a', seq_2a + reverse_complement(seq_2a)])
        amplicon_list.append(['2a_2b', seq_2a + seq_2b])
        amplicon_list.append(['2a_3a', seq_2a + reverse_complement(seq_3a)])
        amplicon_list.append(['2a_3b', seq_2a + seq_3b])

        amplicon_list.append(['2b_2b', reverse_complement(seq_2b) + seq_2b])
        amplicon_list.append(['3a_2b', seq_3a + seq_2b])
        amplicon_list.append(['3b_2b', reverse_complement(seq_3b) + seq_2b])

        amplicon_list.append(['3a_3a', seq_3a + reverse_complement(seq_3a)])
        amplicon_list.append(['3a_3b', seq_3a + seq_3b])

        amplicon_list.append(['3b_3b', reverse_complement(seq_3b) + seq_3b])
        
    write_amplicon(dir_sample, amplicon_info, amplicon_list)



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
# Function to calculate the number of reads and collapsed reads aligned to the reference amplicons
        # this doesn't tell you which type or where
# litterally just gives 2 numbers. Number of amplicons that were aligned and the number after collapse. Doesn't break it down by alignment to any particular amplicon
################################################################################
def analyze_alignments_all_amplicons(dir_sample, amplicon_info, min_MAPQ, min_AS):
    N7 = amplicon_info['index_I1']
    N5 = amplicon_info['index_I2']

    exp_dir = create_filename(dir_sample, N7, N5, 'mainfolder')

    file_UMI = create_filename(dir_sample, N7, N5, 'umifastqgz')
    UMI_dict = create_barcode_dict(file_UMI)
    
    ## THIS was added to get an idea of how many unique reads there are.
    a=list(UMI_dict.values())
    lst2 = [item[0] for item in a]
    print('number of reads before mispriming cleanup',len(lst2))
    print('number of unique umis before mispriming cleanup',len(set(lst2)))
    
    results_folder = os.path.join(exp_dir, 'results')
    if not os.path.exists(results_folder):
        os.mkdir(results_folder)

    results_file = create_filename(dir_sample, N7, N5, 'results_all_amplicons')

    file_sorted_bam_amplicons = create_filename(dir_sample, N7, N5, 'sorted_bam_amplicons')

    bam_in_alignment_file = pysam.AlignmentFile(file_sorted_bam_amplicons, 'rb')

    bam_in_all = bam_in_alignment_file.fetch()

    names_list_amplicons = []
    UMI_list_amplicons = []

    for read in bam_in_all:
        if read.has_tag('AS'):
            read_AS = read.get_tag('AS')
        # We test first if the read is unmapped, otherwise read_AS would be undefined
        if not read.is_unmapped and read.mapping_quality >= min_MAPQ \
                and read_AS >= min_AS and not read.is_secondary:
            UMI_list_amplicons.append(UMI_dict[read.query_name][0])
            names_list_amplicons.append(read.query_name)

    all_amplicons_total_reads = len(set(names_list_amplicons))
    all_amplicons_total_reads_collapsed = len(set(UMI_list_amplicons))

    results_df = pd.DataFrame({'all_amplicons_total_reads': [all_amplicons_total_reads],
                               'all_amplicons_total_reads_collapsed': [all_amplicons_total_reads_collapsed]
                               },
                              columns=['all_amplicons_total_reads',
                                       'all_amplicons_total_reads_collapsed'])

    results_df.to_excel(results_file)

    return results_df
    
# Functions from Uditas to analyze the alignments. These support analyze_alignment()
# This analysis looks at the indel distribution at the different cut sites.

# the follwoing 2 functions- parse_indels() and get_intersection() are used in the 3rd function find_indels().
# parse_indels(), get_intersection(), find_indels() are unchanged from uditas


# find_indels() will be used in analyze_alignments() (along with other functions)


################################################################################
# parse_indels
#
# Author: David Kelly
#
# Parse the CIGAR tuples for the positions and lengths of insertions and
# deletions implied by the aligned read.
#
# Input
#  aligned_read:  pysam AlignedSegment object
#
# Output
#  indels:        list of (position, indel_length) tuples. indel_lengths are
#                  positive for insertions, negative for deletions.
################################################################################
def parse_indels(aligned_read):
    indels = []

    if not aligned_read.is_unmapped and not aligned_read.is_secondary:  # We only look at primary alignments
        ref_i = aligned_read.reference_start
        for operation, length in aligned_read.cigartuples:
            if operation == 0:
                ref_i += length
            elif operation == 1:
                indels.append((ref_i, length))
            elif operation == 2:
                indels.append((ref_i, -length))
                ref_i += length
            else:
                print >> sys.stderr, 'Unrecognized CIGAR operation for %s' % aligned_read.query_name

    return indels


################################################################################
# Function to get the number of intersecting bases between two intervals
################################################################################
def get_intersection(region1_begin, region1_end, region2_begin, region2_end):
    list1 = range(int(region1_begin) + 1, int(region1_end) + 1)
    list2 = range(int(region2_begin) + 1, int(region2_end) + 1)
    return len(set(list1).intersection(list2))


################################################################################
# find_indels
#
# Input
#  bam_file:               alignments file to process
#
################################################################################
def find_indels(bam_file, strand, region_chr, region_start, region_end, UMI_dict, min_MAPQ, min_AS):
    bam_in_alignment_file = pysam.AlignmentFile(bam_file, 'rb')

    # We get the reads that overlap the window in which we make the counts
    # fetch will get reads with any overlap with the window
    # For UDiTaS, care must be take to ensure that the read covers the whole window, some short reads may cover just
    # one side of the window, depending on the direction of the UDiTaS primer
    bam_in = bam_in_alignment_file.fetch(region_chr, region_start, region_end)

    names_list = []
    position_list = []
    indel_list = []
    UMI_list = []
    for read in bam_in:
        # We add here a check to make sure the read came from the primer, crossed the cut and covered the whole window
        if read.has_tag('AS'):
            read_AS = read.get_tag('AS')
        # We test first if the read is unmapped, otherwise read_AS would be undefined
        if not read.is_unmapped and (((read.reference_start < region_start) and
                 (read.reference_end > region_end)) and
                    read.mapping_quality >= min_MAPQ and
                    read_AS >= min_AS and  not read.is_secondary):
            read_indels = parse_indels(read)
            # if no indels found, write 0
            if len(read_indels) == 0:
                read_indels.append(('-', 0))

            # print indels to table
            for pos, indel in read_indels:
                names_list.append(read.query_name)
                if pos == '-':
                    position_list.append(-1)
                else:
                    position_list.append(int(pos))

                indel_list.append(indel)
                UMI_list.append(UMI_dict[read.query_name][0])

    df = pd.DataFrame({'read_name': names_list,
                       'position': position_list,
                       'indel': indel_list,
                       'UMI': UMI_list})

    df['position_end'] = df.position + np.abs(df.indel)

    overlap = [get_intersection(df.loc[index]['position'], df.loc[index]['position_end'], region_start, region_end)
               for index in range(df.shape[0])]

    position_filter = np.array(overlap) > 0

    deletion_filter = position_filter & np.array(df.indel < 0)

    insertion_filter = position_filter & np.array(df.indel > 0)

    total_reads_in_region = len(set(df['read_name']))
    total_collapsed_reads_in_region = len(set(df['UMI']))

    total_indels = len(set(df.loc[position_filter]['read_name']))
    total_collapsed_indels = len(set(df.loc[position_filter]['UMI']))

    total_deletions = len(set(df.loc[deletion_filter]['read_name']))
    total_collapsed_deletions = len(set(df.loc[deletion_filter]['UMI']))

    total_insertions = len(set(df.loc[insertion_filter]['read_name']))
    total_collapsed_insertions = len(set(df.loc[insertion_filter]['UMI']))

    return [total_reads_in_region, total_indels, total_deletions, total_insertions,
            total_collapsed_reads_in_region,
            total_collapsed_indels,
            total_collapsed_deletions,
            total_collapsed_insertions]

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


################################################################################
# Function to analyse fragment sizes
################################################################################
def analyze_fragment_sizes(dir_sample, amplicon_info, min_MAPQ):
    N7 = amplicon_info['index_I1']
    N5 = amplicon_info['index_I2']
    file_sorted_bam_amplicons = create_filename(dir_sample, N7, N5, 'sorted_bam_amplicons')

    # Section to plot
    file_figure_classes = (create_filename(dir_sample, N7, N5, 'results_amplicons') + '_fragment_sizes_classes.pdf')
    file_figure = (create_filename(dir_sample, N7, N5, 'results_amplicons') + '_fragment_sizes.pdf')

    bam_in = pysam.AlignmentFile(file_sorted_bam_amplicons, "rb")

    iter_bam_in = bam_in.fetch()

    segments_bam_in = create_segments(iter_bam_in, bam_in, min_MAPQ)

    df = pd.DataFrame(segments_bam_in, columns=['type', 'begin', 'end'])
    df['length'] = df['end'] - df['begin']

    if df.shape[0] > 0:
        median_size = np.median(df['length'])

        file_fragments = (create_filename(dir_sample, N7, N5, 'results_amplicons') + '_fragment_sizes.xlsx')
        df.to_excel(file_fragments, index=False)

        df2 = df.pivot(columns='type', values='length')

        up_limit = 700
        fs = 20
        plt.rcParams["figure.figsize"] = [20, 10]

        # if df2.shape
        # plot with categories separated by colors
        df2.plot.hist(stacked=True, bins=np.arange(0, up_limit, 20))
        plt.xlim(0, up_limit)
        plt.xlabel('Fragment Size (bp)', fontsize=fs)
        plt.ylabel('Counts', fontsize=fs)
        plt.xticks(fontsize=fs)
        plt.yticks(fontsize=fs)
        if len(plt.gca().get_legend_handles_labels()[0]) > 0:  # To prevent warning from legend
            plt.legend(fontsize=12)
        pylab.savefig(file_figure_classes, bbox_inches='tight')
        plt.close(plt.gcf())

        # plot with all categories with the same color
        plt.hist(df['length'], np.arange(0, up_limit, 20))
        plt.xlim(0, up_limit)
        plt.xlabel('Fragment Size (bp)', fontsize=fs)
        plt.ylabel('Counts', fontsize=fs)
        plt.xticks(fontsize=fs)
        plt.yticks(fontsize=fs)
        if len(plt.gca().get_legend_handles_labels()[0]) > 0:  # To prevent warning from legend
            plt.legend(fontsize=12)
        pylab.savefig(file_figure, bbox_inches='tight')
        plt.close(plt.gcf())
    else:
        median_size = 0

    bam_in.close()

    return median_size


# Functions from Uditas to analyze the alignments. These support analyze_alignment()

# get_cut_in_reference_amplicon_df() is used in znalyze_alignment()
# This had to be altered to get the cuts for REPLACE




#######################################
# Function to get cut list depending on amplicon type
# updated to include HDR
#
#
# Inputs are reaction_type and fasta record
# This function is needed to catch the cases when there are zero (control) or two cuts in the reference amplicon,
# eg when we have two cuts in the same chromosome
#
# The cut positions in the output depend on how the amplicons were constructed in create_amplicon. A change in
# that function must be matched in this function
#
# output: list with pairs of prefixes and cut positions [['cut1', pos1], 'cut2', pos2]]
#######################################
def get_cut_in_reference_amplicon_df(amplicon_info, reaction_type, record, strand, window_size, amplicon_window_around_cut):
    cut_in_reference_amplicon_df = pd.DataFrame(columns=['cut_type', 'cut_position'])
    if reaction_type == 'control':
        if strand == '+':
            # For the control sample (just primer, no cuts) we will look at counts and indels as if there was a cut
            # in position window_size
            cut_in_reference_amplicon_df.loc[0] = ['', window_size + 1]
        else:  # For the - strand we count reads at the end of the amplicon
            cut_in_reference_amplicon_df.loc[0] = ['', amplicon_window_around_cut - window_size - 2]
    elif reaction_type in ['single_cut', 'double_cut_different_chromosomes', 'triple_cut_different_chromosomes']:
        # We need to get cut1 in case its coordinates are smaller than amplicon_window_around_cut
        if amplicon_info['strand_guide_1'] == '+':
            # sp or sa for the moment only
            cut1 = amplicon_info['end_guide_1'] - 3
        elif amplicon_info['strand_guide_1'] == '-':
            cut1 = amplicon_info['start_guide_1'] + 3
        else:
            raise StrandError('strand_guide_1 can only have as values + or -')

        if cut1 < amplicon_window_around_cut:
            cut_site = cut1
        else:
            cut_site = amplicon_window_around_cut

        # all amplicons have a single cut at position amplicon_window_around_cut
        cut_in_reference_amplicon_df.loc[0] = ['cut1', cut_site]
    elif reaction_type in ['double_cut_same_chromosome', 'double_cut_same_chromosome_and_HDR']:
        # Case two guides on the same chromosome
        if amplicon_info['strand_guide_1'] == '+':
            # sp or sa for the moment only
            cut1 = amplicon_info['end_guide_1'] - 3
        elif amplicon_info['strand_guide_1'] == '-':
            cut1 = amplicon_info['start_guide_1'] + 3
        else:
            raise StrandError('strand_guide_1 can only have as values + or -')

        if amplicon_info['strand_guide_2'] == '+':
            cut2 = amplicon_info['end_guide_2'] - 3
        elif amplicon_info['strand_guide_2'] == '-':
            cut2 = amplicon_info['start_guide_2'] + 3
        else:
            raise StrandError('strand_guide_2 can only have as values + or -')

        # We switch the coordinates of cut1 and cut2 if the guides are provided so that cut2 < cut1
        # cut1 it will always be smaller than cut2 and in the results cut1 will be the cut site with
        # smaller genomic coordinate
        # cut1 and cut2 also flipped in create_amplicon

        if cut2 < cut1:
            (cut1, cut2) = (cut2, cut1)
        if cut1 < amplicon_window_around_cut:
            cut_site = cut1
        else:
            cut_site = amplicon_window_around_cut

        # In this case some amplicons (wt, large_inversion) have two cuts, the rest have one
        if record.name in ['wt', 'large_inversion']:
            cut1_cut2_length = cut2 - cut1
            cut_in_reference_amplicon_df.loc[0] = ['cut1', cut_site]
            cut_in_reference_amplicon_df.loc[1] = ['cut2', cut_site + cut1_cut2_length]
        else:
            cut_in_reference_amplicon_df.loc[0] = ['cut1', cut_site]
    elif reaction_type == 'replace':
        # Case for replace targeting
        if amplicon_info['strand_guide_1'] == '+':
            # sp or sa for the moment only
            cut1 = amplicon_info['end_guide_1'] - 3
        elif amplicon_info['strand_guide_1'] == '-':
            cut1 = amplicon_info['start_guide_1'] + 3
        else:
            raise StrandError('strand_guide_1 can only have as values + or -')

        if amplicon_info['strand_guide_2'] == '+':
            cut2 = amplicon_info['end_guide_2'] - 3
        elif amplicon_info['strand_guide_2'] == '-':
            cut2 = amplicon_info['start_guide_2'] + 3
        else:
            raise StrandError('strand_guide_2 can only have as values + or -')

        # We switch the coordinates of cut1 and cut2 if the guides are provided so that cut2 < cut1
        # cut1 it will always be smaller than cut2 and in the results cut1 will be the cut site with
        # smaller genomic coordinate
        # cut1 and cut2 also flipped in create_amplicon

        if cut2 < cut1:
            (cut1, cut2) = (cut2, cut1)
        if cut1 < amplicon_window_around_cut:
            cut_site = cut1
        else:
            cut_site = amplicon_window_around_cut
        
        replace_length = int(len(amplicon_info['Replace_Donor']))
        cut1_cut2_length = cut2 - cut1
        
        # In this case some amplicons (wt, large_inversion) have two cuts, the rest have one
        if record.name in ['wt', 'large_inversion']:
            cut_in_reference_amplicon_df.loc[0] = ['cut1', cut_site]
            cut_in_reference_amplicon_df.loc[1] = ['cut2', cut_site + cut1_cut2_length]            
        elif record.name in ['replace_fwd', 'replace_rev']:
            cut_in_reference_amplicon_df.loc[0] = ['cut1', cut_site]
            cut_in_reference_amplicon_df.loc[1] = ['cut2', cut_site + replace_length]
        #these are for the hiti alignments 
        elif record.name in ['hiti_5prime_replace_exon_fwd', 'hiti_5prime_replace_exon_rev']:
            cut_in_reference_amplicon_df.loc[0] = ['cut1', replace_length]
        elif record.name in ['hiti_3prime_replace_exon_fwd', 'hiti_3prime_replace_exon_rev']:
            cut_in_reference_amplicon_df.loc[0] = ['cut1', cut1_cut2_length]  
        #concatomers
        elif record.name in ['doner_tail_tail', 'doner_head_tail', 'doner_head_head']:
            cut_in_reference_amplicon_df.loc[0] = ['cut1', replace_length]
        else:
            cut_in_reference_amplicon_df.loc[0] = ['cut1', cut_site]
    else:
        raise ReactionTypeError('Reaction type not yet supported by current version of UDiTaS')

    return cut_in_reference_amplicon_df




################################################################################
#  Function to analyze indels and structural rearrangements from aligned reads to amplicons
################################################################################
def analyze_alignments(dir_sample, amplicon_info, window_size, amplicon_window_around_cut, min_MAPQ, min_AS):

    reaction_type = get_reaction_type(amplicon_info)

    # UDiTaS primer strand
    strand = amplicon_info['strand']

    N7 = amplicon_info['index_I1']
    N5 = amplicon_info['index_I2']
    exp_dir = create_filename(dir_sample, N7, N5, 'mainfolder')
    results_folder = os.path.join(exp_dir, 'results')
    if not os.path.exists(results_folder):
        os.mkdir(results_folder)

    results_file = (create_filename(dir_sample, N7, N5, 'results_amplicons') + '_results_amplicon_window_'
                    + str(window_size) + '.xlsx')

    file_sorted_bam_amplicons = create_filename(dir_sample, N7, N5, 'sorted_bam_amplicons')

    file_UMI = create_filename(dir_sample, N7, N5, 'umifastqgz')

    UMI_dict = create_barcode_dict(file_UMI)

    filename_amplicons_fa = os.path.join(exp_dir, 'amplicons', 'amplicons.fa')

    results_df = pd.DataFrame(index=[0])

    # We get the reference amplicon list
    with open(filename_amplicons_fa, "rU") as handle:
        records = list(SeqIO.parse(handle, "fasta"))
    
    #first we get fine the cut sites for each fasta amplicon
    for record in records:
        cut_df = get_cut_in_reference_amplicon_df(amplicon_info, reaction_type, record, strand, window_size,
                                                  amplicon_window_around_cut)

        # We go over cut1 and cut2 and when using cut1 we check if we are in control sample
        for i in cut_df.index:
            cut = cut_df.loc[i, 'cut_type']
            cut_position = cut_df.loc[i, 'cut_position']
            region_chr = record.name

            region_start = cut_position - window_size
            region_end = cut_position + window_size + 1

            results = find_indels(file_sorted_bam_amplicons, strand, region_chr, region_start, region_end, UMI_dict,
                                  min_MAPQ, min_AS)
            # This is to catch the control case, no cut there
            if len(cut) > 0:
                prefix = region_chr + '_' + cut
            else:
                prefix = region_chr

            results_df[prefix + '_total_reads'] = [results[0]]
            results_df[prefix + '_total_indels'] = [results[1]]
            results_df[prefix + '_total_deletions'] = [results[2]]
            results_df[prefix + '_total_insertions'] = [results[3]]
            results_df[prefix + '_total_reads_collapsed'] = [results[4]]
            results_df[prefix + '_total_indels_collapsed'] = [results[5]]
            results_df[prefix + '_total_deletions_collapsed'] = [results[6]]
            results_df[prefix + '_total_insertions_collapsed'] = [results[7]]

    median_size = analyze_fragment_sizes(dir_sample, amplicon_info, min_MAPQ)
    results_df['median_fragment_size'] = [median_size]

    results_df.to_excel(results_file)
    return results_df

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
#       this function is for preparing LAM files made by Misha for  These 'UMIs' are not real. They are just 
#       put together to allow the uditas to work. But there is no UMI in the LAM setup these are random seqs
#
#       The adapter sequence is in read 2 : GACTATAGGGCACGCGTGG
#       both read 1 and 2 have 5nt of random at the begining NNNNN
#
####################################

def umi_joining(dir_sample, amplicon_info):
    length_of_umis_on_end = 5
    adapter_seq = 'GACTATAGGGCACGCGTGG'
    adapt_len = len(adapter_seq)
    
    #define the input files and output files
    N7 = amplicon_info['index_I1']
    N5 = amplicon_info['index_I2']

    #input files
    r1_fastq = create_filename(dir_sample, N7, N5, 'R1fastqgz')
    r2_fastq = create_filename(dir_sample, N7, N5, 'R2fastqgz')
    
    
    #output files not gz compressed #NEED TO FIX!
    file_R1_LAM = create_filename(dir_sample, N7, N5, 'R1fastq_LAM')
    file_R2_LAM = create_filename(dir_sample, N7, N5, 'R2fastq_LAM')
    rUMI_fastq = create_filename(dir_sample, N7, N5, 'umifastq')
    
    
    #open/create all these files
    ref_file_R1_LAM = open(file_R1_LAM, "w")
    ref_file_R2_LAM = open(file_R2_LAM, "w")
    ref_rUMI_fastq = open(rUMI_fastq, "w")
    
    file_count = 0
    
    UMI_list = []
    
    # We open r1,r2 files and make our new files
    with open_fastq_or_gz(r1_fastq) as r1_file, open_fastq_or_gz(r2_fastq) as r2_file:
        # Add counters for all reads

        r1_r2 = itertools.izip(r1_file, r2_file)
        
        
        for header_r1, header_r2 in r1_r2:
    
            seq_r1_UMI, seq_r2_UMI_adapter = r1_r2.next()

            r1_r2.next()

            qual_r1_UMI, qual_r2_UMI_adapter = r1_r2.next()

            seq_r1_UMI, seq_r2_UMI_adapter = seq_r1_UMI.rstrip(), seq_r2_UMI_adapter.rstrip()

            qual_r1_UMI, qual_r2_UMI_adapter = qual_r1_UMI.rstrip(), qual_r2_UMI_adapter.rstrip()
            
            #choose the UMI sequences and qualities. Combin read 1 and read 2
            
            umi1 = seq_r1_UMI[:length_of_umis_on_end]
            umi2 = seq_r2_UMI_adapter[:length_of_umis_on_end]
            umi = umi1 + umi2
            
            UMI_list.append(umi)
            
            umi_qual1 = qual_r1_UMI[:length_of_umis_on_end]
            umi_qual2 = qual_r2_UMI_adapter[:length_of_umis_on_end]
            umi_qual = umi_qual1 + umi_qual2
            
            # take read 1 and 2 and make without the UMI/adapter sequences
            
            seq_r1 = seq_r1_UMI[length_of_umis_on_end:]
            seq_r1_qual = qual_r1_UMI[length_of_umis_on_end:]
            
            length_UMI_adapter = length_of_umis_on_end + adapt_len
            seq_r2 = seq_r2_UMI_adapter[length_UMI_adapter:]
            seq_r2_qual = seq_r2_UMI_adapter[length_UMI_adapter:]
            
            #now print r1, r2, and the umi to fastq files
            
            #read1
            print("\n".join([header_r1.rstrip(), seq_r1.rstrip(), "+", seq_r1_qual.rstrip()]), file=ref_file_R1_LAM)
            #read2
            print("\n".join([header_r1.rstrip(), seq_r2.rstrip(), "+", seq_r2_qual.rstrip()]), file=ref_file_R2_LAM)
            #umis
            print("\n".join([header_r1.rstrip(), umi.rstrip(), "+", umi_qual.rstrip()]), file=ref_rUMI_fastq)
            
            file_count += 1
    
    print('number of reads:', file_count)
    unique_umis = len(set(UMI_list))
    print('number of unique umis', unique_umis)


######################
#
#
#       this is different than my other correct priming in that it looks at the other read.
#
########## this reads a different input file than Tn5!!! #####

def correct_priming_lam(dir_sample, amplicon_info, primer_seq_plus_downstream):
    
    #defined the sequence as input (primer and 12 nt downstream). Normally 32-37nt sequence
    n_max_mismatches = 3
    length_primer_down = len(primer_seq_plus_downstream)
    
    files_out = list()
    n_file = 0
    files_out_dict = dict()

    for file_selected in files_out:
        files_out_dict[os.path.basename(file_selected)] = n_file
        n_file += 1
    
    
    #define the input files and output files
    N7 = amplicon_info['index_I1']
    N5 = amplicon_info['index_I2']
    
    #input files Now altered to input the lam treated files.
    
    r1_fastq = create_filename(dir_sample, N7, N5, 'R1fastq_LAM')
    r2_fastq = create_filename(dir_sample, N7, N5, 'R2fastq_LAM')
    
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

        reads_in_experiment_list_count = 0

        reads_good_prime = 0
        read_misprimed = 0
        read_perfect_match = 0
        read_mismatch_but_good = 0
 
        r1_r2 = itertools.izip(r1_file, r2_file)

        for header_r1, header_r2 in r1_r2:
   
            seq_r1, seq_r2 = r1_r2.next()

            r1_r2.next()

            qual_r1, qual_r2 = r1_r2.next()

            seq_r1, seq_r2 = seq_r1.rstrip(), seq_r2.rstrip()

            qual_r1, qual_r2 = qual_r1.rstrip(), qual_r2.rstrip()

            #We mask with N any bases with scores below or equal to , (11, default in mask)
            
            #don't need seq1
            #seq_r1_use = mask(seq_r1, qual_r1)

            #this is modified for and LAM has the primer on read1
            seq_r1_use = mask(seq_r1[:length_primer_down], qual_r1[:length_primer_down])

    
            # change to 1 for reads with perfect match or match within hamming distance decided above
            is_good_read = 0

            if seq_r1_use == primer_seq_plus_downstream:
                # perfect match case
                is_good_read = 1
                read_perfect_match += 1
            else:
                # We look for barcodes with up to two mismatches, default in select_barcode
                h_d = hamm_dist(seq_r1_use, primer_seq_plus_downstream)
                if h_d <= n_max_mismatches:
                    # match after checking within hamming distance
                    is_good_read = 1
                    read_mismatch_but_good +=1

            reads_in_experiment_list_count += 1

                        
            #Print good reads 
            if is_good_read:
                reads_good_prime += 1
                
                # We test whether the read has on of the combination of indices from our experiment list
                # If not save in a separate file
                r1f = ref_file_R1_corrprime
                r2f = ref_file_R2_corrprime

                print("\n".join([header_r1.rstrip(), seq_r1.rstrip(), "+", qual_r1.rstrip()]), file=r1f)
                print("\n".join([header_r2.rstrip(), seq_r2.rstrip(), "+", qual_r2.rstrip()]), file=r2f)
                    
                

            else:
                # We print reads with mispriming
                print("\n".join([header_r1.rstrip(), seq_r1.rstrip(), "+", qual_r1.rstrip()]), file=ref_file_out_r1_not_corrprime)
                print("\n".join([header_r2.rstrip(), seq_r2.rstrip(), "+", qual_r2.rstrip()]), file=ref_file_out_r2_not_corrprime)
                
                read_misprimed += 1
                
        print('total_read_number', reads_in_experiment_list_count)
        print('reads_with_good_priming', reads_good_prime)
        print('reads_with_bad_priming', read_misprimed)
        print('perfect matched',read_perfect_match)

        print('good but mismatched',read_mismatch_but_good)

 
# functions to align trimmed-to-break reads to the genome


############################
#
#  #input is the correctly primed part from pipeline1
# This is used to trim the fastq files. Originally it was just to trim off the adapter sequences.
# However in this case we use it to trim off the sequences up to the break

# for Tn5 use the reverse compliment and read1. Read 1 starts distal to the site
# for LAM use read2

# ##########################
def trim_fastq_to_break(dir_sample, amplicon_info, seq_primer_to_breaksite, lam_or_tn5='tn5', length=35):

    
    seqRC = reverse_complement(seq_primer_to_breaksite)
    length = str(length) 
    print('minimum length after trimming is:', length)
    # We first check if the experiment had any guides
    N7 = amplicon_info['index_I1']
    N5 = amplicon_info['index_I2']

    # input files
    #for Tn5 we only use read 1 as it is the sequence distal to the gene specific primer
    file_R1 = create_filename(dir_sample, N7, N5, 'R1fastq_CorrPrime')
    file_R2 = create_filename(dir_sample, N7, N5, 'R2fastq_CorrPrime')
    
    if lam_or_tn5 == 'tn5':
        input_file = file_R1
        print('tn5 analysis')
    elif lam_or_tn5 == 'lam':
        input_file = file_R2
        print('lam analysis')
    else:
        print('input has to be "tn5" or "lam" ')
    
     # output files -trimmed to break
    file_cutadapt_R1 = create_filename(dir_sample, N7, N5, 'R1trimmed2break')
    file_cutadapt_R2 = create_filename(dir_sample, N7, N5, 'R2trimmed2break')
    file_cutadapt_report = create_filename(dir_sample, N7, N5, 'trimmed_report2break')
    
    if lam_or_tn5 == 'tn5':
        output_file = file_cutadapt_R1
    elif lam_or_tn5 == 'lam':
        output_file = file_cutadapt_R2
    else:
        print('input has to be "tn5" or "lam" ')
        
    if not os.path.exists(os.path.dirname(file_cutadapt_R1)):
        os.mkdir(os.path.dirname(file_cutadapt_R1))

    # remove adapters with cutadapt Also perhaps increase it to 0.33
    #original uditas peramiter had an error -e 0.33 (but was cutting of random stuff too much)
        #also I extended this to 35 as the 25 alignment was giving a lot of high scoring alignments so I wanted to make sure that they align to a unique genome spot
    cutadapt_command = ['cutadapt',
                        '-m', length,
                        '-e', '0.20',
                        '-a', seqRC,
                        '-o', output_file,
                        input_file]

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
# This quantifies the output of the local or global alignments. It searches for where in the genome
# they alinged and then tallies this, and for the umis
# choosing lam_or_tn5 just changes the output file to contain collaped UMIs or not (LAM doesn't have these)
#
#  input: 
# dir_sample: directory of the file.
# global_or_local: choosing the bam file made by the global or local alignment
# tn5_or_lam is just for organiznig het output plot. Lam has no UMI collapsing so the table just removes the umi info
# targeting_vector_name: this need to be the exact string that was put into the indexed bowtie2 so that it can be found in the alignment file
# window size: this is the region around the cut site imported from the csv file. This region is considered the "target region" 
#                and everything else is the translocation
# min_MAPq is the minimum map quality in the alignment files
# min_AS this is the minimal alignment socre. Not sure what are good numbers but Uditas said -180 for global alignments
#
#
############################


def quantify_pipeline2_alignments(dir_sample, global_or_local, amplicon_info, lam_or_tn5 = 'tn5',
                                  targeting_vector_name = 'pE038_MC',
                                   window_size=10000, min_MAPQ = 25, min_AS=-180):
    #function beginning
    # define samples
    N7 = amplicon_info['index_I1']
    N5 = amplicon_info['index_I2']
    targeted_chromosome = amplicon_info['chr_guide_1']
    
    if global_or_local == 'global':
        bam_file = create_filename(dir_sample, N7, N5, 'break_trimmed_sorted_bam_genome_global')
        results_file = create_filename(dir_sample, N7, N5, 'results_pipeline2_global')
    elif global_or_local == 'local':
        bam_file = create_filename(dir_sample, N7, N5, 'break_trimmed_sorted_bam_genome_local')
        results_file = create_filename(dir_sample, N7, N5, 'results_pipeline2_local')
    
    
    #make the dictionary
    file_UMI = create_filename(dir_sample, N7, N5, 'umifastqgz')
    UMI_dict = create_barcode_dict(file_UMI)

    #make the regions of the chromosome to be anlayzed
    list_of_not_targeted_chromosomes = []
    for i in range(1,23):
        chromsome = ('chr'+str(i))
        if chromsome != targeted_chromosome:
            list_of_not_targeted_chromosomes.append(chromsome)
        elif chromsome == targeted_chromosome:
            #our region of interest
            start_targeting_region = int(amplicon_info['start_guide_1']) - window_size
            end_targeting_region = int(amplicon_info['start_guide_1']) + window_size
    list_of_not_targeted_chromosomes = list_of_not_targeted_chromosomes + ['chrX', 'chrY', 'chrM']


    #prepare counts for measing
    all_reads_unfiltered = []
    all_reads_unfiltered_umis = []
    alignments_target_site_names = []
    alignments_target_site_umis = []
    alignments_targeting_vector_names = []
    alignments_targeting_vector_umis = []
    alignments_on_other_chromosome_regions_names = []
    alignments_on_other_chromosome_regions_umis = []
    alignments_on_other_chromosome_regions_names = []
    
    
    #open up the bam file
    bam_alignment_file = pysam.AlignmentFile(bam_file, 'rb')
    bam_in = bam_alignment_file.fetch()

    #test each read
    for read in bam_in:
        all_reads_unfiltered.append(read.query_name)
        all_reads_unfiltered_umis.append(UMI_dict[read.query_name][0])
        if read.has_tag('AS'):
            read_AS = read.get_tag('AS')
        #check read quality is good
        if read.mapping_quality >= min_MAPQ and read_AS >= min_AS and  not read.is_secondary:
            # should use this as a filtered bam file for plotting
            # count the alignments to the vector
            if read.reference_name == targeting_vector_name:           
                alignments_targeting_vector_names.append(read.query_name)
                alignments_targeting_vector_umis.append(UMI_dict[read.query_name][0])

            elif read.reference_name in list_of_not_targeted_chromosomes:            
                alignments_on_other_chromosome_regions_names.append(read.query_name)
                alignments_on_other_chromosome_regions_umis.append(UMI_dict[read.query_name][0])

            #if alignment matches between start_bin/end_bin then it is the targeted region
            elif read.reference_name == targeted_chromosome:

                if start_targeting_region < read.pos and read.pos < end_targeting_region:
                    alignments_target_site_names.append(read.query_name)
                    alignments_target_site_umis.append(UMI_dict[read.query_name][0])

                #otherwise we add it to the list of translocations if its farther then our targeting region
                else:
                    alignments_on_other_chromosome_regions_names.append(read.query_name)
                    alignments_on_other_chromosome_regions_umis.append(UMI_dict[read.query_name][0])


    ##### calculating the results  ####
    count_all = len(all_reads_unfiltered)
    count_all_umis = len(set(all_reads_unfiltered_umis))
    # targeting site        
    count_target = len(alignments_target_site_names)
    count_target_umis = len(set(alignments_target_site_umis))

    # vector alignments
    count_target_vector = len(alignments_targeting_vector_names)
    count_target_vector_umis = len(set(alignments_targeting_vector_umis))

    # translocations                               
    count_translocations = len(alignments_on_other_chromosome_regions_names)
    count_translocations_umis = len(set(alignments_on_other_chromosome_regions_umis))

    #### totals ###
    total_align = float(count_target + count_target_vector + count_translocations)
    total_align_collapsed = float(count_target_umis + count_target_vector_umis + count_translocations_umis)

    # % of total alignments
    per_target_site = (count_target/total_align)*100
    per_target_vector = count_target_vector/total_align*100
    per_translocations = count_translocations/total_align*100

    # % of collapsed alignments
    per_target_site_collapsed = count_target_umis/total_align_collapsed*100
    per_target_vector_collapsed = count_target_vector_umis/total_align_collapsed*100
    per_translocations_collapsed = count_translocations_umis/total_align_collapsed*100

        

    results_df = pd.DataFrame({'alignments_target_site_all': [count_target],
                               'collapsed_alignments_target_site': [count_target_umis],
                               'alignments_target_vector_all': [count_target_vector],
                               'collapsed_alignments_target_vector': [count_target_vector_umis],
                               'alignments_translocations_all': [count_translocations],
                               'collapsed_alignments_translocations': [count_translocations_umis],
                               'total aligned reads filtered': [total_align],
                               'total collpased reads filtered': [total_align_collapsed],
                               '%_target_site_reads': [per_target_site],
                               '%_target_vector_reads': [per_target_vector],
                               '%_translocations': [per_translocations],
                               '%_target_site_reads_collapsed': [per_target_site_collapsed],
                               '%_target_vector_reads_collapsed': [per_target_vector_collapsed],
                               '%_translocations_collapsed': [per_translocations_collapsed],
                               'all_reads_unfiltered': [count_all],
                               'all_reads_unfiltered_collapsed': [count_all_umis]
                               }, 
                                  columns = ['%_target_site_reads',
                                             '%_target_vector_reads',
                                             '%_translocations', 
                                             '%_target_site_reads_collapsed', 
                                             '%_target_vector_reads_collapsed', 
                                             '%_translocations_collapsed',
                                             'alignments_target_site_all',
                                             'alignments_target_vector_all', 
                                             'alignments_translocations_all',
                                             'collapsed_alignments_target_site',
                                             'collapsed_alignments_target_vector',
                                             'collapsed_alignments_translocations', 
                                             'total aligned reads filtered', 
                                             'total collpased reads filtered',
                                             'all_reads_unfiltered',
                                             'all_reads_unfiltered_collapsed'])
    
    if lam_or_tn5 == 'lam':
        results_df = results_df[['%_target_site_reads',
                                             '%_target_vector_reads',
                                             '%_translocations',
                                             'alignments_target_site_all',
                                             'alignments_target_vector_all', 
                                             'alignments_translocations_all',
                                             'total aligned reads filtered',
                                             'all_reads_unfiltered']]
    
    results_df.to_excel(results_file)

    return results_df


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


###########################################################
#
#               'CORRECT-PRIMING' FILTER UPDATED
#
#   checks for priming on Read1, exports out good sequences a. 
#
#   
#   mismatch = number of mismatches for hamming distance when checing primer or primer_and_downstream
#   exports out a data frame of information about the proccessing as well as the good samples fastq
#   trimmedR1R2 = False if input is raw data, True if you trimed fastq with trimming_R1_R2() and input has been trimmed
#   removePrimerPlusDownstream = True to remove the primer_and_downstream seq used to check file. 
#                                This is best for guideseq or removing integrated exogenous seq.
#                                False keeps the entire sequence intact when exporting.
#   exportMismatch = Fale don't make a mismatch file
#############
def correct_priming2(dir_sample, amplicon_info, primer_seq, primer_seq_plus_downstream, mismatch, trimmed_R1R2, removePrimerPlusDownstream, exportMismatch):
       
    
    length_primer_down = len(primer_seq_plus_downstream)
    length_primer = len(primer_seq)
    n_max_mismatches = mismatch

    
    #define the input files and output files
    N7 = amplicon_info['index_I1']
    N5 = amplicon_info['index_I2']
    name = amplicon_info['description']
    
    #input files
    if trimmed_R1R2 == 0:
        r1_fastq = create_filename(dir_sample, N7, N5, 'R1fastqgz')
        r2_fastq = create_filename(dir_sample, N7, N5, 'R2fastqgz')
    elif trimmed_R1R2 == 1:
        r1_fastq = create_filename(dir_sample, N7, N5, 'R1fastqgz_LAM')
        r2_fastq = create_filename(dir_sample, N7, N5, 'R2fastqgz_LAM')
    
    #output files not gz compressed
    file_R1_corrprime = create_filename(dir_sample, N7, N5, 'R1fastq_CorrPrime')
    file_R2_corrprime = create_filename(dir_sample, N7, N5, 'R2fastq_CorrPrime')

    #mismatched files
    file_out_r1_not_corrprime = os.path.join(dir_sample, N7 + '_' + N5, 'fastq_files', 'mis_prime_R1.fastq')
    file_out_r2_not_corrprime = os.path.join(dir_sample, N7 + '_' + N5, 'fastq_files', 'mis_prime_R2.fastq')
       
    #open/create all these files
    ref_file_R1_corrprime = open(file_R1_corrprime, "w")
    ref_file_R2_corrprime = open(file_R2_corrprime, "w")
    
    if exportMismatch == 1:
        ref_file_out_r1_not_corrprime = open(file_out_r1_not_corrprime, "w")
        ref_file_out_r2_not_corrprime = open(file_out_r2_not_corrprime, "w")
    
    
    # We open r1,r2 files and distribute reads
    with open_fastq_or_gz(r1_fastq) as r1_file, open_fastq_or_gz(r2_fastq) as r2_file:
        # Add counters for all reads

        total_reads = 0
        primer_read = 0
        primer_read_plusdownstream = 0
    
 
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
            
            #export sequence for Read1, either trim off primer for guideseq
            if removePrimerPlusDownstream == 1:
                seq_r1_export = seq_r1[length_primer_down:]
                qual_r1_export = qual_r1[length_primer_down:]
            elif removePrimerPlusDownstream == 0:
                seq_r1_export = seq_r1
                qual_r1_export = qual_r1

            #this is modified for miniseq UMI->Index2
            #seq_r2_use = mask(seq_r2, qual_r1)

            # change to 1 for reads with perfect match or match within hamming distance decided above
            is_good_read = 0
            
            #tests if there is a match for the primer and the nucletides downstream
            if seq_r1_primer_and_downstream == primer_seq_plus_downstream:
                # perfect match case
                is_good_read = 1
            else:
                # We look for barcodes with up to two mismatches, default in select_barcode
                #h_d = hamm_dist(seq_r1_primer_and_downstream, primer_seq_plus_downstream)
                v_d = damerauLevenshtein(seq_r1_primer_and_downstream, primer_seq_plus_downstream, similarity=False)
                #if h_d <= n_max_mismatches:
                if v_d <= n_max_mismatches:
                    # match after checking within hamming distance
                    is_good_read = 1
            
            #tests if there is a match for the primer (could have or not half downstream)
            if seq_r1_primer_size == primer_seq:
                primer_read += 1
            else:
                #h_d = hamm_dist(seq_r1_primer_size, primer_seq)
                v_d = damerauLevenshtein(seq_r1_primer_size, primer_seq, similarity=False)
                #if h_d <= n_max_mismatches:
                if v_d <= n_max_mismatches:
                    primer_read += 1
            
            #Print good reads 
            if is_good_read:
                primer_read_plusdownstream += 1
                
                # We test whether the read has on of the combination of indices from our experiment list
                # If not save in a separate file
                r1f = ref_file_R1_corrprime
                r2f = ref_file_R2_corrprime

                print("\n".join([header_r1.rstrip(), seq_r1_export.rstrip(), "+", qual_r1_export.rstrip()]), file=r1f)
                print("\n".join([header_r2.rstrip(), seq_r2.rstrip(), "+", qual_r2.rstrip()]), file=r2f)
                    
        
            elif exportMismatch == 1:
                # We print reads with mispriming
                print("\n".join([header_r1.rstrip(), seq_r1.rstrip(), "+", qual_r1.rstrip()]), file=ref_file_out_r1_not_corrprime)
                print("\n".join([header_r2.rstrip(), seq_r2.rstrip(), "+", qual_r2.rstrip()]), file=ref_file_out_r2_not_corrprime)
                
        
    results_file = create_filename(dir_sample, N7, N5, 'localpriming_data')
    
    
    results_df = pd.DataFrame({'sample_name': [name],
                               'i7': [N7],
                               'i5': [N5],
                               'total_reads': [total_reads],
                               'reads_with_good_priming': [primer_read_plusdownstream],
                               'reads_with_gene_specific_primer_misprimed': [primer_read - primer_read_plusdownstream],
                               'reads_with_indexing_primer_mispriming': [total_reads - primer_read], #indexing primer i
                               },
                                  columns=['sample_name',
                                           'i7',
                                           'i5',
                                           'total_reads',
                                           'reads_with_good_priming',
                                           'reads_with_gene_specific_primer_misprimed',
                                           'reads_with_indexing_primer_mispriming'])
    
    results_df.to_excel(results_file)
    return results_df



############################   TRIM OFF ADAPTERS IN SHORT AMPLICONS ##############
#
#      This uses trim adapter software to remove the ends of the amplicons if they were too short
#           and then the sequence was running onto Illumina adapters or onto the LAM ligated adapters
#
#      In this I use the known primer sequneces to cut the adapter (the NNNNN go after in LAM so I trim before)
#      And I use the adapter sequence of the LAM for the 3' end of Read1. 
#                    4)I could use the illumina adapters also as an alternative strategy


def trim_short_fastq(dir_sample, amplicon_info, direction5primer, direction3primer, adapter):

    
    direction = amplicon_info['Direction']

    # Read1  This is the gene specific primering binding
    if direction == 5:
        gene_specific_primer = direction5primer         # this is the gene specific primer for the 5 side
    elif direction == 3:
        gene_specific_primer = direction3primer         # gene specific primer for the 3 side
    rev_gene_specific_primer = reverse_complement(gene_specific_primer)
    
    # LAM adapter
    adapter_rev = reverse_complement(adapter)
    
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
                        '-a', adapter_rev,
                        '-A', rev_gene_specific_primer,
                        '-o', file_cutadapt_R1, '-p', file_cutadapt_R2,
                        file_R1, file_R2]

    handle_cutadapt_report = open(file_cutadapt_report, 'wb')
    subprocess.call(cutadapt_command, stdout=handle_cutadapt_report)
    handle_cutadapt_report.close()
    




    
################################################
#           Function to write reference AAV sequence
# 
#     if using HDR arms in the AAV then you might want to remove this sequence 
#    
################################################
def create_AAV_reference(dir_sample, amplicon_info):
    N7 = amplicon_info['index_I1']
    N5 = amplicon_info['index_I2']
    exp_dir = create_filename(dir_sample, N7, N5, 'mainfolder')
    amplicon_folder = os.path.join(exp_dir, 'amplicons')
    if not os.path.exists(amplicon_folder):
        os.mkdir(amplicon_folder)


    # If several plasmids were given, separate them using ';' in sample_info.csv. Here we remove the ';' so that
    # we just concatenate the sequences
    if len(amplicon_info['AAVseq']) < 5:
        print('there is no AAV sequence')
        
        return

    filename = os.path.join(exp_dir, amplicon_folder, 'AAV.fa')
    file_handle = open(filename, "w")    
    pl_seq = amplicon_info['AAVseq'].replace(';', '')
    seq1 = Seq(pl_seq, IUPAC.unambiguous_dna)
    record1 = SeqRecord(seq1, 'AAV', description='')
    SeqIO.write(record1, file_handle, 'fasta')

    file_handle.close()
    # Create index file
    initial_dir = os.getcwd()
    os.chdir(amplicon_folder)
    index_err_file = os.path.join(amplicon_folder, 'index_AAV.err')
    index_out_file = os.path.join(amplicon_folder, 'index_AAV.out')

    index_err_fh = open(index_err_file, 'wb')
    index_out_fh = open(index_out_file, 'wb')
    subprocess.call(['bowtie2-build',
                     filename, 'AAV'], stderr=index_err_fh, stdout=index_out_fh)
    os.chdir(initial_dir)
    index_err_fh.close()
    index_out_fh.close() 
    

def extract_reads_plasmid(dir_sample, amplicon_info):

    N7 = amplicon_info['index_I1']
    N5 = amplicon_info['index_I2']
    
    #input sorted read
    file_sorted_bam_plasmid_local = create_filename(dir_sample, N7, N5, 'sorted_bam_plasmid_local')
    
    #UNMAPPED FILES
    #file to output unmapped reads to
    file_unmapped_bam_plasmid = create_filename(dir_sample, N7, N5, 'unmapped_bam_plasmid_local')
    
    #file for sorting unmapped reads
    file_qsorted_unmapped_bam_plasmid = create_filename(dir_sample, N7, N5, 'qsorted_unmapped_bam_plasmid_local')

    #files for exporting unmapped fastq
    file_R1_unmapped = create_filename(dir_sample, N7, N5, 'unmapped_plasmid_R1fastq')
    file_R2_unmapped = create_filename(dir_sample, N7, N5, 'unmapped_plasmid_R2fastq')

    if not os.path.exists(os.path.dirname(file_R1_unmapped)):
        os.mkdir(os.path.dirname(file_R1_unmapped))
    
    # MAPPED PLASMID FILES
    #file to output Mapped plasmid reads to
    file_mapped_bam_plasmid = create_filename(dir_sample, N7, N5, 'mapped_bam_plasmid_local')
    
    #file for sorting unmapped reads
    file_qsorted_mapped_bam_plasmid = create_filename(dir_sample, N7, N5, 'qsorted_mapped_bam_plasmid_local')

    #files for exporting unmapped fastq
    file_R1_mapped = create_filename(dir_sample, N7, N5, 'mapped_plasmid_R1fastq')
    file_R2_mapped = create_filename(dir_sample, N7, N5, 'mapped_plasmid_R2fastq')
    
    if not os.path.exists(os.path.dirname(file_R1_mapped)):
        os.mkdir(os.path.dirname(file_R1_mapped))
    
    #EXTRACTING UNMAPPED
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
    
   
    #EXTRACTING MAPPED
    #         -G function means skip if it contains all flags, 0xC =12 which is 3 and 4 unmapped and unmaped pair
    extract_mapped_bam_command = ['samtools', 'view', '-b', '-G', '0xC', file_sorted_bam_plasmid_local, '-o',
                                    file_mapped_bam_plasmid]

    subprocess.call(extract_mapped_bam_command)

    qsort_mapped_bam_command = ['samtools', 'sort', '-n', file_mapped_bam_plasmid, '-o',
                                  file_qsorted_mapped_bam_plasmid]

    subprocess.call(qsort_mapped_bam_command)

    bamtofastq_command = ['bedtools', 'bamtofastq', '-i', file_qsorted_mapped_bam_plasmid,
                          '-fq', file_R1_mapped, '-fq2', file_R2_mapped]

    file_err2 = file_R1_mapped[:-9] + '_err.txt'
    handle_file_err = open(file_err2, 'wb')

    subprocess.call(bamtofastq_command, stderr=handle_file_err)

    for fo in [file_R1_mapped, file_R2_mapped]:
        with open(fo) as f_in, gzip.open(fo + '.gz', 'wb') as f_out:
            f_out.writelines(f_in)
        os.remove(fo)

    

########## Aligns the extracted fastq files from plasmid align ##########
#        
#
#
        
def align_AAV_local(dir_sample, amplicon_info, ncpu=4):

    # We first check if the experiment had any guides
    N7 = amplicon_info['index_I1']
    N5 = amplicon_info['index_I2']
    # exp_dir = create_filename(dir_sample, N7, N5, 'mainfolder')
    
    #need to open up 'mapped_plasmid_R1fastqgz' because it was a zipped file
    file_cutadapt_R1 = create_filename(dir_sample, N7, N5, 'mapped_plasmid_R1fastqgz')
    file_cutadapt_R2 = create_filename(dir_sample, N7, N5, 'mapped_plasmid_R2fastqgz')
        
    file_sam_AAV_local = create_filename(dir_sample, N7, N5, 'sam_AAV_local')
    file_sam_report_AAV_local = create_filename(dir_sample, N7, N5, 'sam_report_AAV_local')

    if not os.path.exists(os.path.dirname(file_sam_AAV_local)):
        os.mkdir(os.path.dirname(file_sam_AAV_local))
             
    file_bam_AAV_local = create_filename(dir_sample, N7, N5, 'bam_AAV_local')
    file_sorted_bam_AAV_local = create_filename(dir_sample, N7, N5, 'sorted_bam_AAV_local')
    # file_sorted_bai_genome_local = create_filename(dir_sample, N7, N5, 'sorted_bai_genome_local')

    if not os.path.exists(os.path.dirname(file_bam_AAV_local)):
        os.mkdir(os.path.dirname(file_bam_AAV_local))        
          
    print('made filenames')    
        
    # local alignment to the genome with bowtie2
    initial_dir = os.getcwd()

    folder_amplicons = create_filename(dir_sample, N7, N5, 'amplicons')

    os.chdir(folder_amplicons)

    # --local alignment, -X space -k number hit
    bowtie2_command = ['bowtie2', '--local', '-p', str(ncpu),
                       '-X', '4000', '-k', '2', '-x', 'AAV',
                             '-1', file_cutadapt_R1, '-2', file_cutadapt_R2,
                             '-S', file_sam_AAV_local]

    handle_sam_report_genome_local = open(file_sam_report_AAV_local, 'wb')

    subprocess.call(bowtie2_command, stderr=handle_sam_report_genome_local)

    handle_sam_report_genome_local.close()

    # convert sam to bam
    sam_to_bam_AAV_local_command = ['samtools', 'view', '-Sb', file_sam_AAV_local]

    handle_file_bam_AAV_local = open(file_bam_AAV_local, 'wb')

    subprocess.call(sam_to_bam_AAV_local_command, stdout=handle_file_bam_AAV_local)

    # sort bam files
    sort_bam_AAV_local_command = ['samtools', 'sort', file_bam_AAV_local, '-o', file_sorted_bam_AAV_local]

    subprocess.call(sort_bam_AAV_local_command)

    # Create bam index files
    create_bam_AAV_local_index_command = ['samtools', 'index', file_sorted_bam_AAV_local]
    subprocess.call(create_bam_AAV_local_index_command)

    # Clean up
    #os.remove(file_sam_AAV_local)
    #os.remove(file_bam_AAV_local)

    os.chdir(initial_dir)
    

#analyze plasmid alignments

# analysis = 'plasmid' is the analysis of alignment to the plasmid
# min mapQ should be >1

def analyze_local_alignments(dir_sample, amplicon_info, min_MAPQ, analysis):
    N7 = amplicon_info['index_I1']
    N5 = amplicon_info['index_I2']
    name = amplicon_info['name']
    
    exp_dir = create_filename(dir_sample, N7, N5, 'mainfolder')

    results_folder = os.path.join(exp_dir, 'results')
    if not os.path.exists(results_folder):
        os.mkdir(results_folder)

    results_file = create_filename(dir_sample, N7, N5, 'results_plasmid')

    if analysis == 'plasmid':
        file_sorted_bam_plasmid_local = create_filename(dir_sample, N7, N5, 'sorted_bam_plasmid_local')

        bam_in_alignment_file = pysam.AlignmentFile(file_sorted_bam_plasmid_local, 'rb')
        bam_in = bam_in_alignment_file.fetch()

        names_list_local_plasmid = []
        
        for read in bam_in:
            if read.mapping_quality >= min_MAPQ and not read.is_unmapped and not read.is_secondary:
                names_list_local_plasmid.append(read.query_name)
            
        total_reads_plasmid = len(set(names_list_local_plasmid))

        results_df = pd.DataFrame({'sample_name': [name],
                               'i7': [N7],
                               'i5': [N5],
                               'number_plasmid_alignments': [total_reads_plasmid],
                                },
                                  columns=['sample_name',
                                           'i7',
                                           'i5',
                                           'number_plasmid_alignments'])

    results_df.to_excel(results_file)

    return results_df    




################################################################################
#  Function to analyze indels and structural rearrangements from aligned reads to amplicons

# In htis I removed the UMI analysis. Also I added calculating number of aligned reads 
#    and I calculate HDR by any amplicon that mapps well to the HDR sequence (there is not clear cut site)
################################################################################

def analyze_alignments_LAM(dir_sample, amplicon_info, window_size, amplicon_window_around_cut, min_MAPQ, min_AS):

    reaction_type = get_reaction_type(amplicon_info)

    # UDiTaS primer strand
    strand = amplicon_info['strand']

    N7 = amplicon_info['index_I1']
    N5 = amplicon_info['index_I2']
    name = amplicon_info['name']
    description = amplicon_info['description']

    
    exp_dir = create_filename(dir_sample, N7, N5, 'mainfolder')
    results_folder = os.path.join(exp_dir, 'results')
    if not os.path.exists(results_folder):
        os.mkdir(results_folder)

    results_file = (create_filename(dir_sample, N7, N5, 'results_amplicons') + '_results_amplicon_window_'
                    + str(window_size) + '.xlsx')

    file_sorted_bam_amplicons = create_filename(dir_sample, N7, N5, 'sorted_bam_amplicons')


    filename_amplicons_fa = os.path.join(exp_dir, 'amplicons', 'amplicons.fa')

    results_df = pd.DataFrame(index=[0])
    HDRcount = 0
    mappedreads = 0
    # We get the reference amplicon list
    with open(filename_amplicons_fa, "rU") as handle:
        records = list(SeqIO.parse(handle, "fasta"))
    
    #Add name and index of the read
    results_df['sample_name'] = [name]
    results_df['i7'] = [N7]
    results_df['i5'] = [N5]
    results_df['sample_descript'] = [description]

    #count all the good alignments
    bam_in_alignment_file = pysam.AlignmentFile(file_sorted_bam_amplicons, 'rb')
    bam_in = bam_in_alignment_file.fetch()
    for read in bam_in:
        min_MAPQ = 4

        if read.has_tag('AS'):
            read_AS = read.get_tag('AS')
                   # We test first if the read is unmapped, otherwise read_AS would be undefined
        if not read.is_unmapped and read.mapping_quality >= min_MAPQ and read_AS >= min_AS and not read.is_secondary:
            mappedreads += 1
    print('Allreads is : ', mappedreads)
    results_df['Total_mapped' + '_reads'] = [mappedreads]
    
    
    #first we get fine the cut sites for each fasta amplicon
    for record in records:
        cut_df = get_cut_in_reference_amplicon_df(amplicon_info, reaction_type, record, strand, window_size,
                                                  amplicon_window_around_cut)

        if record.name == 'HDR_knockin':

            min_MAPQ = 4
            
            bam_in_alignment_file = pysam.AlignmentFile(file_sorted_bam_amplicons, 'rb')
            bam_in = bam_in_alignment_file.fetch('HDR_knockin')

            

            for read in bam_in:
                if read.has_tag('AS'):
                        read_AS = read.get_tag('AS')
                    # We test first if the read is unmapped, otherwise read_AS would be undefined
                if not read.is_unmapped and read.mapping_quality >= min_MAPQ and read_AS >= min_AS and not read.is_secondary:
                    HDRcount += 1
            print('HDRcount is : ', HDRcount)
            results_df['HDR' + '_total_reads'] = [HDRcount]

    # We get the reads that overlap the window in which we make the counts
    # fetch will get reads with any overlap with the window
    # For UDiTaS, care must be take to ensure that the read covers the whole window, some short reads may cover just
    # one side of the window, depending on the direction of the UDiTaS primer

        else:
            # We go over cut1 and cut2 and when using cut1 we check if we are in control sample
            for i in cut_df.index:
                cut = cut_df.loc[i, 'cut_type']
                cut_position = cut_df.loc[i, 'cut_position']
                region_chr = record.name

                region_start = cut_position - window_size
                region_end = cut_position + window_size + 1

                results = find_indelsLAM(file_sorted_bam_amplicons, strand, region_chr, region_start, region_end, min_MAPQ, min_AS)
                # This is to catch the control case, no cut there
                if len(cut) > 0:
                    prefix = region_chr + '_' + cut
                else:
                    prefix = region_chr

                results_df[prefix + '_total_reads'] = [results[0]]
                results_df[prefix + '_total_indels'] = [results[1]]
                results_df[prefix + '_total_deletions'] = [results[2]]
                results_df[prefix + '_total_insertions'] = [results[3]]

    median_size = analyze_fragment_sizes(dir_sample, amplicon_info, min_MAPQ)
    results_df['median_fragment_size'] = [median_size]

    results_df.to_excel(results_file)
    return results_df








################################################################################
# find_indels  modified to remove the UMI portion
#
# Input
#  bam_file:               alignments file to process
#
################################################################################
def find_indelsLAM(bam_file, strand, region_chr, region_start, region_end, min_MAPQ, min_AS):
    bam_in_alignment_file = pysam.AlignmentFile(bam_file, 'rb')

    # We get the reads that overlap the window in which we make the counts
    # fetch will get reads with any overlap with the window
    # For UDiTaS, care must be take to ensure that the read covers the whole window, some short reads may cover just
    # one side of the window, depending on the direction of the UDiTaS primer
    bam_in = bam_in_alignment_file.fetch(region_chr, region_start, region_end)

    names_list = []
    position_list = []
    indel_list = []
    for read in bam_in:
        # We add here a check to make sure the read came from the primer, crossed the cut and covered the whole window
        if read.has_tag('AS'):
            read_AS = read.get_tag('AS')
        # We test first if the read is unmapped, otherwise read_AS would be undefined
        if not read.is_unmapped and (((read.reference_start < region_start) and
                 (read.reference_end > region_end)) and
                    read.mapping_quality >= min_MAPQ and
                    read_AS >= min_AS and  not read.is_secondary):
            read_indels = parse_indels(read)
            # if no indels found, write 0
            if len(read_indels) == 0:
                read_indels.append(('-', 0))

            # print indels to table
            for pos, indel in read_indels:
                names_list.append(read.query_name)
                if pos == '-':
                    position_list.append(-1)
                else:
                    position_list.append(int(pos))

                indel_list.append(indel)
               

    df = pd.DataFrame({'read_name': names_list,
                       'position': position_list,
                       'indel': indel_list})

    df['position_end'] = df.position + np.abs(df.indel)

    overlap = [get_intersection(df.loc[index]['position'], df.loc[index]['position_end'], region_start, region_end)
               for index in range(df.shape[0])]

    position_filter = np.array(overlap) > 0

    deletion_filter = position_filter & np.array(df.indel < 0)

    insertion_filter = position_filter & np.array(df.indel > 0)

    total_reads_in_region = len(set(df['read_name']))
    total_indels = len(set(df.loc[position_filter]['read_name']))
    total_deletions = len(set(df.loc[deletion_filter]['read_name']))
    total_insertions = len(set(df.loc[insertion_filter]['read_name']))
    
    
    
    return [total_reads_in_region, 
            total_indels, 
            total_deletions, 
            total_insertions]


############################
#
# Aligns reads globally to amplicon. "end-to-end" as the default function in bowtie2
# Input: directory to be analyzed
#        amplicon_info, slice of sample_info.csv for the sample being processed
#        file_genome_2bit, 2bit file with the reference genome being used
#
#        paired = True or False
# ##########################
def align_ampliconLAM(dir_sample, amplicon_info, check_plasmid_insertions, paired, ncpu=4):

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
    
    if paired == True:
        bowtie2_command = ['bowtie2', '-p', str(ncpu), '--very-sensitive',
                           '-X', '5000', '-k', '2', '-x', 'amplicons',
                           '-1', file_R1, '-2', file_R2,
                           '-S', file_sam_amplicons]
    
    else:
        bowtie2_command = ['bowtie2', '--very-sensitive', '-p', str(ncpu),
                   '-X', '5000', '-k', '2', '-x', 'amplicons',
                   '-U', file_R1, '-S', file_sam_amplicons]


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
    #os.remove(file_sam_amplicons)
    os.remove(file_bam_amplicons)

    # Create bam index files
    create_bam_amplicons_index_command = ['samtools', 'index', file_sorted_bam_amplicons]
    subprocess.call(create_bam_amplicons_index_command)

    os.chdir(initial_dir)
 
    
    
    
    
    ################ Trimming off x nucleotides from the 5' ends #####
#   inputs files, trims off a defined amount for 5' ends of read 1 or two and then exports compressed file.
#
# R1trim, <int> the length to trim  read1
# R2trim, <int> the length to trim read2
# pipeline, '1' if trimming fastq that are normal R1fastqgz. 
#           '2' if for pipeline 2 and fastq extracted from plasmid reads
#      
#  for pipeline 1 this is used to chop off the 5' ends of sequences, uses the unsorted fastq files
#  for pipeline 2 this is to trim up to the break and uses the primer checked, non-plasmid aligned fastq

def trimming_R1_R2(dir_sample, amplicon_info, R1trim, R2trim, pipeline):
    
    
    #define the input files and output files
    N7 = amplicon_info['index_I1']
    N5 = amplicon_info['index_I2']
    
    #input/output files
    if pipeline == 1:
        #inputs
        r1_fastq = create_filename(dir_sample, N7, N5, 'R1fastqgz')
        r2_fastq = create_filename(dir_sample, N7, N5, 'R2fastqgz')
        
         #output files not gz compressed
        r1_clipped_fastq = create_filename(dir_sample, N7, N5, 'R1fastq_LAM')
        r2_clipped_fastq = create_filename(dir_sample, N7, N5, 'R2fastq_LAM')
    
    elif pipeline == 2:
        print('makingpipeline2 files')
        r1_fastq = create_filename(dir_sample, N7, N5, 'unmapped_plasmid_R1fastqgz')
        r2_fastq = create_filename(dir_sample, N7, N5, 'unmapped_plasmid_R2fastqgz')
        
        #output files not gz compressed
        r1_clipped_fastq = create_filename(dir_sample, N7, N5, 'breaktrimmed_R1fastq')
        r2_clipped_fastq = create_filename(dir_sample, N7, N5, 'breaktrimmed_R2fastq')
   

    
    #open/create all these files
    ref_file_r1_fasta_good = open(r1_clipped_fastq, "w")
    ref_file_r2_fasta_good = open(r2_clipped_fastq, "w")

      
    
    # We open r1,r2 files and distribute reads
    with open_fastq_or_gz(r1_fastq) as r1_file, open_fastq_or_gz(r2_fastq) as r2_file:

        r1_r2 = itertools.izip(r1_file, r2_file)

        for header_r1, header_r2 in r1_r2:
           
            seq_r1, seq_r2 = r1_r2.next()

            r1_r2.next()

            qual_r1, qual_r2 = r1_r2.next()
            seq_r1, seq_r2 = seq_r1.rstrip(), seq_r2.rstrip()
            qual_r1, qual_r2 = qual_r1.rstrip(), qual_r2.rstrip()

            seq_r1_use = seq_r1[R1trim:]
            seq_q1_use = qual_r1[R1trim:]
            seq_r2_use = seq_r2[R2trim:]
            seq_q2_use = seq_r2[R2trim:]
            
            r1f = ref_file_r1_fasta_good
            r2f = ref_file_r2_fasta_good

            print("\n".join([header_r1.rstrip(), seq_r1_use.rstrip(), "+", seq_q1_use.rstrip()]), file=r1f)
            print("\n".join([header_r2.rstrip(), seq_r2_use.rstrip(), "+", seq_q2_use.rstrip()]), file=r2f)
                    
           
    for fo in [r1_clipped_fastq, r2_clipped_fastq]:
        with open(fo) as f_in, gzip.open(fo + '.gz', 'wb') as f_out:
            f_out.writelines(f_in)
        os.remove(fo)
