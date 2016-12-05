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
    parser.add_argument("-SEreads", "--single end reads", help=".fastq file \
    containing single end RNA seq reads", type=str, required=False)
    parser.add_argument("-PEreads", "--paired end reads", help=".fastq file \
    containing paired end RNA seq reads", type=str, required=False)
    parser.add_argument("-clip", "--ILLUMINACLIP", help="file which contains \
    the sequences that need to be trimmed", type=str, required=False)


def run_trimmomatic(fastq_filename, out_filename, clipper_file, LEADING=3, 
                    TRAILING=3, win_size=4, req_qual=30):
    """ Returns the output of the annotation tool 'augustus' to a file
    
    Keyword arguments:
        fastq_filename -- file in fastq format (.fastq) 
        out_filename -- string, filename for the outputfile
        clipper_file -- string, filename with adaptersequences that need to 
        be trimmed of the raw fastq reads
        LEADING -- integer, specifies minimum quality required to keep a base.
        TRAILING -- integer, specifies minimum quality required to keep a base.
        win_size -- integer, specifies the number of bases to average across
        req_qual -- integer, specifies the average quality required
        SLIDINGWINDOW -- <
    Returns:
        out_filename -- .fastq, name of output file in .fastq format
    """
    
    # checks whether the file already exists, if not it runs the tool
    if os.path.exists(fileout_augustus):
        print "The file already exists"
        return fileout_augustus
    else:
        cmd = "augustus --genemodel=%s --species=%s %s > %s"\
               %(genemodel_option, SPECIES, input_fasta, fileout_augustus)
      
    output_check = subprocess.check_output(cmd, shell=True)
    call_check = subprocess.check_call(cmd, shell = True)
    return call_check #must be 0
      


    

if __name__ == "__main__":
    #Get input file names from command line
    get_arguments() 
    # run trimmomatic tool from command line
    run_augustus(yeast_fasta, "complete", "saccharomyces_cerevisiae_S288C",
                 "augustus.gff")