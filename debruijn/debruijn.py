#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 22 11:22:33 2019

@author: mleprohon
"""

import os
import sys
import getopt
import pytest
import pylint
import networkx

def args_check(argv):
    """This function collect parameters from the input command line, and check that every needed
        parameter is provided and correct. Then it returns a list of parameter or call use().
        Parameter:
            argv : a list of the different argument in the input commande line
        Output:
            arg_list : a list of parameter values 
    """
    # Check that every needed argument are provided.
    try:
        opts, args = getopt.getopt(argv,
                                   "i:k:o:",
                                   ["init_file=", "kmer_size=", "config="])
    except getopt.GetoptError:
        sys.exit(2)

    # Initiation of parameters and then assignation.
    config_file = None
    fastq_file = None
    kmer_size = 21

    for opt, arg in opts:                   
        if opt in ("-i", "--init_file"):
            fastq_file = arg
        elif opt in ("-k", "--kmer_sizes"):
            kmer_size = arg
        elif opt in ("-o", "--config"):
            config_file = arg       
        
            
    arg_list = [fastq_file, kmer_size, config_file]

    return arg_list

def read_fastq(fastq_file):
    with open(fastq_file, "r") as file_in:
        for line in file_in:
            yield next(file_in)
            next(file_in)
            next(file_in)
            
def cut_kmer(seq, k):
    for i in range(len(seq) - k + 1):
        yield seq[i:i+k]    

def main():
    parameters = args_check(sys.argv[1:])
    read_list = read_fastq(parameters[0])
    for i in read_list:
        print(i)
        for j in cut_kmer(i[:-1], 3):
            print(j)    

if __name__ == "__main__":
    main()