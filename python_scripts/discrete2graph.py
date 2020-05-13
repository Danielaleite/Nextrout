import os
import pickle as pkl
import decimal
import click
import sys
import subprocess
from decimal import Decimal
import networkx as nx
from shutil import copyfile
from scipy.spatial import ConvexHull, convex_hull_plot_2d
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection
import matplotlib.pyplot as plt
import time
#--------------------------------------------------
import source_sink_generator
import filtering
import terminal_computation
import pre_extraction
import utils
#-------------------------------------------------

"""
Scripts
"""


def graph_filtering_from_dat_files(folder_name, t, graph_type, beta_d, funct, min_, btns_factor_source, btns_factor_sink,
                                   weighting_method, weighting_method_simplification, source_flag, sink_flag,
                                   BP_weights, reduction_flag):
    '''
    This script computes the filtering of a Graph pre-extracted from the outputs of the DMK solver.
    :param folder_name: folder path where the outputs will be stored. It should be written "./runs/folder_name".
    :param t: threshold for the weights of the edges after filtering.
    :param graph_type: 1 (to use edges and vertices of the grid), 2 (to use only edges), 3 (to use elements
    of the grid).
    :param funct: weights to be assigned to the edges. Either 'tdens' or 'flux'.
    :param min_: threshold for the weights of the edges after pre-extraction.
    :param btns_factor_source: threshold for the nodes in the source region (see more in paper).
    :param btns_factor_sink: threshold for the nodes in the sink region (see more in paper).
    :param weighting_method: 'ER', 'AVG'.
    :param weighting_method_simplification: 'ER', 'IBP', 'BPW'.
    :param source_flag: flag used to define the source region in the continuous problem.
    :param sink_flag: flag used to define the sink region in the continuous problem.
    :param BP_weights: 'BPtdens' to use optimal transport density as weights for edges, 'BPflux' to use optimal flux.
    :param reduction_flag: If 'yes', then the filtered graph is reduced using bifurcation_paths().
    :return:
        G_final_simplification: filtered graph (networkx graph).
    '''
    '''
    This script computes the simplification (both 1st level and 2nd level) of a Graph extracted from the
    continuous solution of the DMK equations.
    Inputs:
    - (folder_name,t,graph_type,funct,weighting_method), -- extracted graph,
    - min_, -- threshold (%) for the edges in the simplification,
    - btns_factor, -- upper bound for the betweenness centrality. Nodes with a greater btns centrality than this value, are not taken as terminals
    - weighting_method_simplification, -- ['IBP','ER','BPW']
        - IBP: ignore BP weights, i.e., keep the original ones,
        - ER: effective reweighting,
        - BPW: BP weighting.
    - source_flag,sink_flag, -- used for plotting reasons,
    - BP_weights, -- opt_tdens or opt_flux weights for the simplification
    '''

    # No simplification for the case below:

    if graph_type == '3' and weighting_method == 'ER':
        print('this wms does not apply for graph type 3.')
    else:
        # Defining the content of the weight.dat file

        weight_flag = 'length'

        # Loading the Graph

        t_string = '%.0E' % decimal.Decimal(str(t))
        if os.getcwd().split('/')[-1] == 'simplifications':
            path_ = './' + folder_name[2:] + '/' + funct + '/' + funct + '_t' + t_string + '_graph' + str(
                graph_type) + '_wm' + str(weighting_method) + '.dat'
        else:

            path_ = '../muffe_sparse_optimization/simplifications/' + folder_name[
                                                                      2:] + '/' + funct + '/' + funct + '_t' + t_string + '_graph' + str(
                graph_type) + '_wm' + str(weighting_method) + '.dat'
        with open(path_, 'rb') as file:
            Graph = pkl.load(file)

        # Computing cc

        #components = list(nx.connected_components(Graph))
        terminal_info = [source_flag, sink_flag, btns_factor_source, btns_factor_sink]

        G_final_simplification, newGraph, ncc, possible_terminals_source, possible_terminals_sink, mapping = filtering.filtering(Graph,
                                                                                                                          folder_name + '/' + funct,
                                                                                                                       beta_d,
                                                                                                                          graph_type,
                                                                                                                          weighting_method_simplification,
                                                                                                                          weight_flag,
                                                                                                                          min_,
                                                                                                                          BP_weights,
                                                                                                                          terminal_info)


        # Plotting the 1st simplification

        btns_factor_source_string = '%.0E' % decimal.Decimal(str(btns_factor_source))
        btns_factor_sink_string = '%.0E' % decimal.Decimal(str(btns_factor_sink))
        min_string = '%.0E' % decimal.Decimal(str(min_))
        t_string = '%.0E' % decimal.Decimal(str(t))

        # First plot

        fig, ax = plt.subplots(1, 1, figsize=(25, 25))
        posG = nx.get_node_attributes(Graph, 'pos')
        posGsimpl = nx.get_node_attributes(G_final_simplification,'pos')

        frame = nx.path_graph(5)
        pos_frame = {0: (0.01, 0.01),
                     1: (0.01, .99),
                     2: (.99, .99),
                     3: (.99, 0.01),
                     4: (0.01, 0.01)}
        nx.draw(frame, pos_frame, node_size=1, width=1, edge_color='white', node_color='b', ax=ax)

        nx.draw_networkx(Graph, posG, node_size=100, width=1, with_labels=False, edge_color='blue', alpha=0.4,
                         node_color='blue', ax=ax)
        nx.draw_networkx(G_final_simplification, posGsimpl, node_size=150, width=4, with_labels=False, edge_color='red',
                         node_color='red', ax=ax)

        terminals = nx.Graph()
        original_label_terminals_list = []
        pos_terminals = {}
        node_label = 1
        for i in range(1, ncc + 1):
            for node in list(possible_terminals_source[i].union(possible_terminals_sink[i])):
                original_label_terminals_list.append(mapping[i][int(node)])
                terminals.add_node(node_label)
                x, y = newGraph[i].node[node]['pos']
                pos_terminals[node_label] = (x, y)
                node_label += 1
        print('len', len(terminals.nodes()))

        if source_flag in ['1_rect', '3_rect', '1_obstacles', '3_obstacles', '2_rect', '2_obstacles', 'center']:
            ax = source_sink_generator.source_plot(source_flag, ax)
        else:
            source_sink_generator.source_plot(source_flag)
        if sink_flag in ['1_rect', 'circle']:
            ax = source_sink_generator.sink_plot(sink_flag, ax)
        else:
            source_sink_generator.sink_plot(sink_flag)

        nx.draw_networkx_nodes(terminals, pos_terminals, node_size=200, edge_color='g', node_color='g', ax=ax)

        plt.title('Original graph (blue) / simplification (red)')
        path_fig = '../muffe_sparse_optimization/simplifications/' + folder_name[
                                                                      2:] + '/' + funct + '/BP_orig_simpl_' + funct + '_t' + t_string + '_graph' + str(
            graph_type) + '_min' + min_string + '_btns_so_' + btns_factor_source_string + '_btns_si_' + btns_factor_sink_string + '_wm' + str(
            weighting_method) + '_wms' + str(weighting_method_simplification)
        plt.axis('on')
        fig.savefig(path_fig + '.png')
        plt.show(block=False)
        time.sleep(.5)
        plt.close('all')

        # Second plot

        plt.figure(1, figsize=(25, 20))
        plt.figure(1)

        edges_, weights_ = zip(*nx.get_edge_attributes(G_final_simplification, 'weight').items())

        ec = nx.draw_networkx_edges(G_final_simplification, posGsimpl, edge_color=weights_, width=2.,
                                    edge_cmap=plt.cm.jet)
        nc = nx.draw_networkx_nodes(G_final_simplification, posGsimpl, node_color='white',
                                    with_labels=False, node_size=100, cmap=plt.cm.jet)
        nc.set_edgecolor('k')

        if source_flag in ['1_rect', '3_rect', '1_obstacles', '3_obstacles', '2_rect', '2_obstacles', 'center']:
            ax = source_sink_generator.source_plot(source_flag, ax)
        else:
            source_sink_generator.source_plot(source_flag)
        if sink_flag in ['1_rect', 'circle']:
            ax = source_sink_generator.sink_plot(sink_flag, ax)
        else:
            source_sink_generator.sink_plot(sink_flag)

        nx.draw(frame, pos_frame, node_size=1, width=1, edge_color='white', node_color='white')

        nx.draw_networkx_nodes(terminals, pos_terminals, node_size=100, edge_color='g', node_color='g')
        plt.title('Weighted simplification: original (' + weighting_method + ')/ simpl (' + str(
            weighting_method_simplification) + ')')
        plt.colorbar(ec)
        plt.axis('on')

        path_ = '../muffe_sparse_optimization/simplifications/' + folder_name[
                                                                      2:] + '/' + funct + '/BP_simpl_' + funct + '_t' + t_string + '_graph' + str(
            graph_type) + '_min' + min_string + '_btns_so_' + btns_factor_source_string + '_btns_si_' + btns_factor_sink_string + '_wm' + str(
            weighting_method) + '_wms' + str(weighting_method_simplification)
        plt.savefig(path_ + '.png')
        plt.show(block=False)
        time.sleep(.5)
        plt.close('all')
        with open(path_ + '.dat', 'wb') as file:
            pkl.dump(G_final_simplification, file)

        # Last simplification
        if reduction_flag == 'yes':
            G_final_simplification = relabeling(G_final_simplification, Graph)
            posGsimpl = nx.get_node_attributes(G_final_simplification, 'pos')
            G_final_simplification_reduction = bifurcation_paths(G_final_simplification, original_label_terminals_list)
            deg = nx.degree_centrality(G_final_simplification_reduction)
            posGsimplred = nx.get_node_attributes(G_final_simplification_reduction, 'pos')

            fig, ax = plt.subplots(1, 1, figsize=(25, 25))
            posG = nx.get_node_attributes(Graph, 'pos')
            frame = nx.path_graph(5)
            pos_frame = {0: (0.01, 0.01),
                         1: (0.01, .99),
                         2: (.99, .99),
                         3: (.99, 0.01),
                         4: (0.01, 0.01)}
            nx.draw(frame, pos_frame, node_size=1, width=1, edge_color='white', node_color='b', ax=ax)

            nx.draw_networkx(G_final_simplification, posGsimpl, node_size=100, width=1, with_labels=False,
                             edge_color='blue', alpha=0.4, node_color='blue', ax=ax)
            nx.draw_networkx(G_final_simplification_reduction, posGsimplred, node_size=150, width=4, with_labels=False,
                             edge_color='red', node_color='red', ax=ax)
            '''
            terminals=nx.Graph()
            pos_terminals={}
            node_label=1
            for i in range(1,len(newGraphList)+1):
                for node in list(possible_terminals_source[i].union(possible_terminals_sink[i])):
                    terminals.add_node(node_label)
                    x, y = newGraph[i].node[node]['pos']
                    pos_terminals[node_label]=(x,y)
                    node_label += 1
            print('len',len(terminals.nodes()))

            '''

            if source_flag in ['1_rect', '3_rect', '1_obstacles', '3_obstacles', '2_rect', '2_obstacles', 'center']:
                ax = source_sink_generator.source_plot(source_flag, ax)
            else:
                source_sink_generator.source_plot(source_flag)
            if sink_flag in ['1_rect', 'circle']:
                ax = source_sink_generator.sink_plot(sink_flag, ax)
            else:
                source_sink_generator.sink_plot(sink_flag)

            nx.draw_networkx_nodes(terminals, pos_terminals, node_size=200, edge_color='g', node_color='g', ax=ax)

            plt.title('Original graph (blue) / simplification (red)')
            path_fig = '../muffe_sparse_optimization/simplifications/' + folder_name[
                                                                      2:] + '/' + funct + '/BP_orig_2nd_simspl_' + funct + '_t' + t_string + '_graph' + str(
                graph_type) + '_min' + min_string + '_btns_so_' + btns_factor_source_string + '_btns_si_' + btns_factor_sink_string + '_wm' + str(
                weighting_method) + '_wms' + str(weighting_method_simplification)
            plt.axis('on')
            fig.savefig(path_fig + '.png')
            plt.show(block=False)
            time.sleep(.5)
            plt.close('all')
            with open(path_ + '.dat', 'wb') as file:
                pkl.dump(G_final_simplification, file)

            return G_final_simplification


