#!/usr/bin/env python
"""
Created on Wed Dec 14 11:08:05 2016

@author: kersj001

collecting only gene names iwth complete name
"""

from __future__ import division
from sys import argv

def only_genes(gene_list):
    
    genes = []
    
    # skip the header line
    gene_list.readline()
    
    for line in gene_list:
        line = line.strip()
        line = line.split(',')
        genes.append(line[1])
    
    return genes

def complete_name(genes):
    """ In our output the "T" is missing between the 'CRO_' and the number, 
    which contain in the database genelist
    """
    
    gene_list = []
    
    for gene in genes:
        gene = gene.split('_')
        gene = '_T'.join(gene)
        gene_list.append(gene)
       
    return gene_list
    

def order_genes(genes):
    
    gene_list = sorted(genes)
    return gene_list

   
def make_new_list(gene_list):    
    output_file = open('GeneList.txt', 'w')
    
    for gene in gene_list:
        line = '{}\n'.format(gene)
        output_file.write(line)
        
        
if __name__ == "__main__":
    file_list = argv[1]
    open_file_list = open(file_list, "r")
    gene_list = only_genes(open_file_list)
    complete_gene_list = complete_name(gene_list) 
    sorted_gene_list = order_genes(complete_gene_list) 
    make_new_list(sorted_gene_list)