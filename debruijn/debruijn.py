#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 22 11:22:33 2019

@author: mleprohon
"""
import statistics
import os
import sys
import getopt
import pytest
import pylint
import networkx as nx
import random

def args_check(argv):
    """ This function collect parameters from the input command line, and check that every needed
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
    """ This function reads a fastq file and returns it as an iterativ object
        Parameter:
            fastq_file : a string representing a path to a file
    """
    with open(fastq_file, "r") as file_in:
        for line in file_in:
            yield next(file_in).strip()
            next(file_in)
            next(file_in)
            
def cut_kmer(seq, k):
    """ This function cut a sequence into k-mers and returns them as an 
        iterative object.
        Parameter:
            seq : a string representing a sequence.
            k : an int representing the length of the extracted words
    """
    for i in range(len(seq) - k + 1):
        yield seq[i:i+k]    

def build_kmer_dict(fastq_file, k):
    """ This function read a fastq file and returns a dictionnary of all 
        k-mers in sequences in it with their abundance.
        Parameter:
            fastq_file : a string representing a path to a file
            k : an int representing the length of the extracted words
    """
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
    """ This function build a graph with all kmers (of size k-1)
        weigthed by their abundance.
        Parameters:
            kmer_dict : a dictionnary of k-mers with their abundance
    """
    tree_graph = nx.DiGraph()
    for i in kmer_dict:
        tree_graph.add_edge(i[:-1], i[1:], weight=kmer_dict[i])
    return tree_graph
 
def get_starting_nodes(tree_graph):
    """ This function returns a list of the starting nodes in a graph.
        Paramaters:
            tree_graph : an nx.Graph
    """
    starting_nodes = []
    for node in tree_graph.nodes:
        if len(tree_graph.pred[node]) == 0:
            starting_nodes.append(node)
    return starting_nodes

def get_sink_nodes(tree_graph):
    """ This function returns a list of the sink nodes in a graph.
        Paramaters:
            tree_graph : an nx.Graph
    """
    starting_nodes = []
    for node in tree_graph.nodes:
        if len(tree_graph.succ[node]) == 0:
            starting_nodes.append(node)
    return starting_nodes

def get_contigs(tree_graph, starting_nodes, ending_nodes):
    """ This function returns a list of all contigs from a De Bruyn Graph
        of k-mers by looking at all possible path between starting and sink 
        nodes.
        Parameters:
            tree_graph : an nx.Graph representing a De Bruyn Graph of k-mers.
            starting_nodes : a list of starting node.
            sink_nodes : a list of sink nodes.
    """
    contig_list = []
    for i in starting_nodes:
        for j in ending_nodes:
            path_ite = nx.algorithms.simple_paths.all_simple_paths(
                     tree_graph,
                     i,
                     j
                     )
            for k in path_ite:
                 contig = k[0] 
                 for l in range(1, len(k)):
                     contig = contig + k[l][-1]
                 contig_tuple = (contig, len(contig))
                 contig_list.append(contig_tuple)
    return contig_list

def fill(text, width=80):
    """ Split text with a line return to respect fasta format.
        Parameters:
            text : a string representing a text
            width : an int representing the formatting line width.
    """
    return os.linesep.join(text[i:i+width] for i in range(0, len(text), width))

def save_contigs(contig_list, output_file):
    """ This function save a list of contigs in a file in fasta format.
        Parameters:
            contig_list : a list of contigs.
            output_file : a string representing a path to and output file.
    """
    with open(output_file, "w") as filout:
        for i in range(len(contig_list)):
            filout.write(">contig_{} len={}\n".format(i, contig_list[i][1]))
            filout.write(fill(contig_list[i][0]) + "\n")

def std(val_list):
    """ Call statistics.stdev on a list of values.
        Parameter:
            val_list : a list of values.
    """
    return statistics.stdev(val_list)

def path_average_weight(graph, graph_path):
    """ This fonction compute the average weigth of edge in a graph path.
        Parameters:
            graph : a nx.Graph.
            graph_path : a path in graph.
    """
    weight_list = []
    for u, v, e in graph.subgraph(graph_path).edges(data=True):
        weight_list.append(e['weight'])
    return statistics.mean(weight_list)    


def remove_paths(
        graph,
        graph_paths,
        delete_entry_node=False,
        delete_sink_node=False
        ):
    """ This function removes nodes belonging to a path in graph.
        Parameters:
            graph : a nx.Graph . 
            graph_paths : a list of path in graph.
            delete_entry_node : Boolean, delete starting node of a
            path if True.
            delete_sink_node : Boolean, delete sink node of a
            path if True.
    """
    for i in range(len(graph_paths)):
        node_to_remove = list(graph_paths[i])
        if not delete_entry_node:
            node_to_remove = node_to_remove[1:]
        if not delete_sink_node:
            node_to_remove = node_to_remove[:-1]
        for n in node_to_remove:
            graph.remove_node(n)
        return graph

def select_best_path(
        graph,
        graph_paths,
        graph_path_lengths,
        graph_path_weights,
        delete_entry_node=False,
        delete_sink_node=False,
        ):
    """ This function removes nodes belonging to a path in graph.
        Parameters:
            graph : a nx.Graph . 
            graph_paths : a list of path in graph.
            graph_path_lengths : a list of average length of paths in graph.
            graph_path_weights : a list of average weights of paths in graph.
            delete_entry_node : Boolean, delete starting node of a
            path if True.
            delete_sink_node : Boolean, delete sink node of a
            path if True.
    """
    max_weight = max(graph_path_weights)
    candidate = []
    for i in range(len(graph_path_weights)):
        if graph_path_weights[i] == max_weight:
            candidate.append(i)
    if len(candidate) > 1:
        max_length = 0
        best = []
        for i in candidate:
            if graph_path_lengths[i] >= max_length:
                max_length = graph_path_lengths[i]
                best.append(i)
        if len(best) > 1:
            best_path = random.sample(best, 1)
            best_path = best_path[0]
        else:
            best_path = best[0]
    else:
        best_path = candidate[0]
    bad_paths = graph_paths[:best_path] + graph_paths[best_path+1:]
    return remove_paths(graph, bad_paths, delete_entry_node, delete_sink_node)

def main():
    parameters = args_check(sys.argv[1:])
    kmer_dict = build_kmer_dict(parameters[0], parameters[1]) 
    tree_graph = build_graph(kmer_dict)
    starting_nodes = get_starting_nodes(tree_graph)
    ending_nodes = get_sink_nodes(tree_graph)
    contig_list = get_contigs(tree_graph, starting_nodes, ending_nodes)
    save_contigs(contig_list, parameters[2])
    graph_4 = nx.DiGraph()
    graph_4.add_edges_from([(1, 2), (3, 2), (2, 4), (4, 5), (2, 8), (8, 9),
                                (9, 5), (5, 6), (5, 7)])
    graph_4 = select_best_path(graph_4, [[2, 4, 5], [2, 8, 9, 5]],
                                          [1, 4], [10, 10])
    return 0

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