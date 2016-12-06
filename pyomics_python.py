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
    parser.add_argument("-SEr", nargs="*", "--se_reads", help=".fastq file \
    containing single end RNA seq reads", type=str, required=False)
    parser.add_argument("-PEr", nargs="*", "--pe_reads", help=".fastq file \
    containing paired end RNA seq reads", type=str, required=False)
    parser.add_argument("-clip", "--illumina_clip", help="file containing \
    the sequences that need to be trimmed off the target sequence", 
    type=str, required=False)


def run_trimmomatic(fastq_filename, clipper_file, lead_val=3, 
                    tail_val=3, win_size=4, req_qual=30):
    """ Returns the trimmed reads of single or paired end data to a file
    
    Keyword arguments:
        fastq_filename -- file with  RNA-seq reads in fastq format (.fastq) 
        trimm_outfile -- string, filename for the outputfile
        clipper_file -- string, filename with adaptersequences that need to 
        be trimmed of the raw fastq reads
        LEADING -- integer, specifies minimum quality required to keep a base.
        TRAILING -- integer, specifies minimum quality required to keep a base.
        win_size -- integer, specifies the number of bases to average across
        req_qual -- integer, specifies the average quality required
    Returns:
        trimm_outfile -- .fastq, name of trimmed RNA-seq reads in .fastq format
    """
    
    trimm_outfile = "trimmed_%s.fastq"%(fastq_filename)
    # checks whether the file already exists, if not it runs the tool
    if os.path.exists(trimm_outfile):
        print "The file already exists"
        return trimm_outfile
    else:
        cmdSE = "TrimmomaticSE %s %s ILLUMINACLIP:%s LEADING:%d TRAILING:%d \
                        SLIDINGWINDOW:%d:%d"
                %(fastq_filename, trimm_outfile, clipper_file, lead_val, 
                  trail_val, win_size, req_qual)
        cmdPE = "TrimmomaticPE %s %s ILLUMINACLIP:%s LEADING:%d TRAILING:%d \
                        SLIDINGWINDOW:%d:%d"
                %(fastq_filename, trimm_outfile, clipper_file, lead_val, 
                  trail_val, win_size, req_qual)
      
    output_check = subprocess.check_output(cmd, shell=True)
    call_check = subprocess.check_call(cmd, shell = True)
    return call_check #must be 0
      


    

if __name__ == "__main__":
    #Get input file names from command line
    arguments = get_arguments() 
    # run trimmomatic tool from command line
    run_trimmomatic(arguments.se_reads, arguments.illumina_clip)