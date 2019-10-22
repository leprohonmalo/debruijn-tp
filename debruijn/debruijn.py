#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 22 11:22:33 2019

@author: mleprohon
"""

import os
import sys
import getopt
import math
import pytest
import pylint
import networkx as nx

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
            kmer_size = int(arg)
        elif opt in ("-o", "--config"):
            config_file = arg       
        
            
    arg_list = [fastq_file, kmer_size, config_file]

    return arg_list

def read_fastq(fastq_file):
    with open(fastq_file, "r") as file_in:
        for line in file_in:
            yield next(file_in).strip()
            next(file_in)
            next(file_in)
            
def cut_kmer(seq, k):
    for i in range(len(seq) - k + 1):
        yield seq[i:i+k]    

def build_kmer_dict(fastq_file, k):
    read_list = read_fastq(fastq_file) 
    kmer_dict = {}
    for i in read_list:
        for kmer in cut_kmer(i, k):
            if kmer in kmer_dict:
                kmer_dict[kmer] = kmer_dict[kmer] + 1
            else:
                kmer_dict[kmer] = 1
    return kmer_dict

def build_graph(kmer_dict):
    tree_graph = nx.DiGraph()
    for i in kmer_dict:
        tree_graph.add_edge(i[:-1], i[1:], weight=kmer_dict[i])
    return tree_graph
 
def get_starting_nodes(tree_graph):
    starting_nodes = []
    for node in tree_graph.nodes:
        if len(tree_graph.pred[node]) == 0:
            starting_nodes.append(node)
    return starting_nodes

def get_sink_nodes(tree_graph):
    starting_nodes = []
    for node in tree_graph.nodes:
        if len(tree_graph.succ[node]) == 0:
            starting_nodes.append(node)
    return starting_nodes

def get_contigs(tree_graph, starting_nodes, ending_nodes):
    contig_list = []
    for i in starting_nodes:
        for j in ending_nodes:
             path_ite = nx.algorithms.simple_paths.all_simple_paths(tree_graph,
                                                        i,
                                                        j)
             for k in path_ite:
                 contig = k[0] + k[1:][-1]
                 contig_tuple = (contig, len(contig))
                 contig_list.append(contig_tuple)
    return contig_list

def fill(text, width=80):
    """Split text with a line return to respect fasta format"""
    return os.linesep.join(text[i:i+width] for i in range(0, len(text), width))

def save_contigs(contig_list, output_file):
    with open(output_file, "w") as filout:
        for i in range(len(contig_list)):
            filout.write("contig_{} len = {}\n".format(i, contig_list[i][1]))
            filout.write(fill(contig_list[i][0]) + "\n")

def std(val_list):
    pass

def main():
    output = "my_contigs.fasta"
    parameters = args_check(sys.argv[1:])
    kmer_dict = build_kmer_dict(parameters[0], parameters[1]) 
    tree_graph = build_graph(kmer_dict)
    starting_nodes = get_starting_nodes(tree_graph)
    ending_nodes = get_sink_nodes(tree_graph)
    print("START\n", starting_nodes)
    print("END\n", ending_nodes)
    contig_list = get_contigs(tree_graph, starting_nodes, ending_nodes)
    print(contig_list)
    save_contigs(contig_list, output)
    return 0

    

def path_average_weight():
    pass


def remove_paths():
    pass


def select_best_path():
    pass


def solve_bubble():
    pass


def simplify_bubbles():
    pass


def solve_entry_tips():
    pass


def solve_out_tips():
    pass

if __name__ == "__main__":
    main()