#!/bin/usr/env python
"""
Created on Mon Dec  5 09:15:07 2016
author: team members of PyOmics project group
Pipeline script for doing a Differential Expression analysis
"""

from __future__ import division
from sys import argv
import argparse
import subprocess
import os.path
import re



def get_arguments():
    """ Gets the arguments from the command line

    Keyword arguments:
        no parameters for this function
    """
    parser = argparse.ArgumentParser(description="Pipeline script for DE \
                                                    analysis")
    parser.add_argument("-SEr", "--se_reads", nargs="*", help=".fastq file \
    containing single end RNA seq reads", type=str, required=False)
    parser.add_argument("-PEr", "--pe_reads", nargs="*", help=".fastq file \
    containing paired end RNA seq reads", type=str, required=False)
    parser.add_argument("-clip", "--illumina_clip", help="file containing \
    the sequences that need to be trimmed off the target sequence", 
    type=str, required=False)
    arguments = parser.parse_args()
    return arguments


def run_trimmomatic(fastq_filename, clipper_file, seed_mm=2, palin_th=30,
                    simple_th=10, lead_val=3, trail_val=3, win_size=4, 
                    req_qual=30):
    """ Returns the trimmed reads of single or paired end data to a file
    
    Keyword arguments:
        fastq_filename -- file with  RNA-seq reads in fastq format (.fastq) 
        trimm_outfile -- string, filename for the outputfile
        clipper_file -- string, filename with adaptersequences that need to 
        be trimmed of the raw fastq reads
        seed_mm -- integer, seedMismatches: specifies maximum mismatch count
        palin_th -- integer, palindromeClipThreshold: specifies how accurate
        the match between two adapter ligated reads must be.
        simple_th -- integer, simpleClipThreshold: specifies how accurate the 
        match between any adapter must be against a read
        LEADING -- integer, specifies minimum quality required to keep a base.
        TRAILING -- integer, specifies minimum quality required to keep a base.
        win_size -- integer, windowSize:specifies the number of bases to 
        average across
        req_qual -- integer, requiredQuality: specifies the average quality 
        required
        single -- boolean, 
    Returns:
        trimm_outfile -- .fastq, name of trimmed RNA-seq reads in .fastq format
    """
    # make for loop to go through every dataset
    
    trimm_outfile = "trimmed_%s"%(fastq_filename)
    # checks whether the file already exists, if not it runs the tool
    if os.path.exists(trimm_outfile):
        print "The file already exists"
        return trimm_outfile
    else:
        cmdSE = "TrimmomaticSE %s %s ILLUMINACLIP:%s:%d:%d:%d LEADING:%d \
        TRAILING:%d SLIDINGWINDOW:%d:%d" %(fastq_filename, trimm_outfile, 
                                           clipper_file, seed_mm, palin_th, 
                                           simple_th, lead_val, trail_val, 
                                           win_size, req_qual)
                                           
        output_check = subprocess.check_output(cmdSE, shell=True)
        call_check = subprocess.check_call(cmdSE, shell=True)
        return call_check #must be 0

#==============================================================================
#         cmdPE = "TrimmomaticPE %s %s ILLUMINACLIP:%s:%d:%d:%d LEADING:%d \
#         TRAILING:%d SLIDINGWINDOW:%d:%d" %(fastq_filename, trimm_outfile, 
#                                            clipper_file, seed_mm, palin_th, 
#                                            simple_th, lead_val, trail_val, 
#                                            win_size, req_qual)
#                            
#         output_check = subprocess.check_output(cmdPE, shell=True)
#         call_check = subprocess.check_call(cmdPE, shell=True)
#         return call_check #must be 0
#==============================================================================

def run_hisat2(ref_genome, splicesites, trimmed_input):
    """
    
    """
    hisat_outfile = "%s_mapping.sam"%(ref_genome)
    # checking if file already exists
    if os.path.exists(hisat_outfile):
        print "The file already exists"
        return hisat_outfile
    else:
        cmd = "hisat2 -x %s --known-splicesites-infile %s -U %s \
        --dta-cfufflinks -S %s"%(ref_genome, splicesites, trimmed_input, 
                                 hisat_outfile)
                                 
        output_check = subprocess.check_output(cmd, shell=True)
        call_check = subprocess.check_call(cmd, shell=True)
        return call_check #must be 0
    


if __name__ == "__main__":
    #Get input file names from command line
    arguments = get_arguments() 
    print arguments
    
    # running the trimmomatic tool
    # only runs when se_reads contain something
    if arguments.se_reads:
        for i in range(len(arguments.se_reads)):
            run_trimmomatic(arguments.se_reads[i], arguments.illumina_clip)
        
#==============================================================================
#         
#     if arguments.pe_reads==True:
#         for i in range(len(arguments.pe_reads)):
#             run_trimmomatic(arguments.pe_reads[i], arguments.illumina_clip,
#                             single=False)
#                             
#     if arguments.se_reads==True:
#         for i in range(len(arguments.se_reads)):
#             run_trimmomatic(arguments.se_reads[i], arguments.illumina_clip, 
#                             single=True)
#                             
#                     
#         
#     # run trimmomatic tool from command line
#     run_trimmomatic(arguments.se_reads[0], arguments.illumina_clip)
#==============================================================================
