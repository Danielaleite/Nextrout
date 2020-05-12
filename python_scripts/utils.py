import networkx as nx
import numpy as np
import itertools
from decimal import Decimal
import matplotlib.pyplot as plt
import pickle as pkl
import contextlib
import os
from Getting_sources_and_sinks import *
import time
import decimal
from source_sink_generator import *


""" 
Changes to the current working folder.
"""

'''
@contextlib.contextmanager
def change_into(dir):
    """goes into dir and go back to the original directory afterwards."""
    current_dir = os.getcwd()

    try:
        os.chdir(dir)
        yield
    except Exception as e:
        raise e
    finally:
        os.chdir(current_dir)


def execute(operations_in_dir):
    """
    Executes commands in a directory.
    Keys are directories, values are lists of commands.
    """

    def execute_in_dir(dir, commands):
        with change_into(dir):
            for cmd in commands:
                os.system(cmd)

    for dir, commands in operations_in_dir.items():
        execute_in_dir(dir, commands)
'''

""" ________________________________________________________________________________________________________________
"""


def get_baryc(folder_name):
    """
    This script imports the graph structure. Its takes only the barycenter positions out of the graph_cell; no edges.
    """
    if os.getcwd().split("/")[-1] == "simplifications":
        file1_ = open("../" + folder_name[2:] + "input/graph_cell.dat", "r")
    else:
        file1_ = open(folder_name + "input/graph_cell.dat", "r")

    graph_coord_triang = file1_.readlines()
    file1_.close()
    n_nodes = int(graph_coord_triang[0][:12])
    graph_coordinates = graph_coord_triang[2 : 2 + n_nodes]
    return graph_coord_triang, graph_coordinates, n_nodes


"""
def bar2dict(graph_coordinates, n_nodes):
    
    Saving the positions of the barycenters in a dict(to plot!)
    
    bar_pos = {}
    n = 1
    bar_pos_list = [line[:-1].split(' ') for line in graph_coordinates]
    # print(bar_pos_list)
    for line_ in bar_pos_list:
        # print(line_)
        bar_pos[str(n)] = []
        line_ = list(dict.fromkeys(line_))
        line_.remove('')
        # print(line_)
        if len(line_) < 3:
            line_.pop(1)
            line_.append(line_[0])
        else:
            line_.pop(2)
        for i in line_:
            bar_pos[str(n)].append(float(i))
        n += 1
    for key_ in bar_pos.keys():
        bar_pos[key_] = np.array(bar_pos[key_])
    return bar_pos
    # print(bar_pos.keys())
"""

# utils!
def extracting_weights(folder_name, file):
    """
    Extracting the weights (tdens)
    -file: string containing the weights
    """
    #

    file_ = open(folder_name + "/" + file, "r")
    file_weights = file_.readlines()
    file_.close()

    return file_weights


# utils!
def weight2dict(file_weights, file_name):
    """
    Saving the weights in a dict
    -Each key corresponds to the "number" of the vertex. The number was given by Enrico's file
    """
    # Cleaning the data, extracting the WEIGHTS and storing them in a dict, converting them into floats
    dict_weights_func = {}
    weights = []
    print("FILE_NAME", file_name)
    n = 1
    firstIndex = 0
    if (
        file_name == "output/result/opt_tdens.dat"
        or file_name == "output/result/opt_tdens.dat"
        or file_name == "output/result/opt_nrm_grad_avg.dat"
        or file_name == "output/result/opt_pot.dat"
    ):
        firstIndex = 3
    # elif file_name == "output/result/opt_flux.dat":
    #    firstIndex=0

    weights_ = [line[:-1].split(" ") for line in file_weights[firstIndex:]]
    for line_ in weights_:
        line_ = [elem for elem in line_ if elem != ""]
        """
        if file_name=="output/result/opt_tdens.dat" or file_name=="output/result/opt_nrm_grad_avg.dat"  or file_name== "opt_nrm_grad_avg.dat"  \
        or file_name=="tdens" or file_name=='flux' or file_name=="opt_tdens.dat" :
            line_.remove('')
        """
        if line_[0] != "time":
            if len(line_) > 1:
                num_ = float(line_[1])
                dict_weights_func[int(line_[0])] = num_
            else:
                num_ = float(line_[0])
                dict_weights_func[n] = num_

            weights.append(num_)
        n += 1
    return dict_weights_func, weights


def completing_with_zeros(dict_weights_func, weights, bar_pos):
    """
    Completing with zeros
    """
    list_ = dict_weights_func.keys()
    # print(list_)
    for n in range(1, len(bar_pos.keys()) + 1):
        if n not in list_:
            dict_weights_func[int(str(n))] = 0
    return dict_weights_func


def func_as_weighted_graph(folder_name, file_name):
    """
    Turning a function in a dat file into a weighted graph,
    where the nodes correspond to the barycentres of the triangles
    in the grid and the weights to the value of the function in
    that particular triangle
        Inputs:
        -folder_name: it is the path to the folder that contains the files (string!)
                    ex: '../../../data/Runs_08_07_2019_0_refinements/runs/run_delta_1/'
        -file_name: it is the name of the file. It could be source, sink, tdens0, 
                    opt_tdens,... (string!)
    """
    plotting = "yes"

    """Importing the graph structure (we take only the bareycenters)."""
    graph_coord_triang, graph_coordinates, n_nodes = get_baryc(folder_name)
    """Saving the positions of the barycenters in a dict(to plot!)"""
    bar_pos = bar2dict(graph_coordinates)
    """Extracting the weights (source)"""
    file_weights = extracting_weights(folder_name, file_name)
    """Saving the source in a dict"""
    dict_weights_func, weights = weight2dict(file_weights, file_name)
    """Completing with zeros"""
    dict_weights_func = completing_with_zeros(dict_weights_func, weights, bar_pos)
    """Defining the graph with the attributes (pos and source)"""
    X_func, pos = dict2graph(bar_pos, dict_weights_func)
    if plotting == "yes":
        """Defining the size of the nodes"""
        size_func = node_size(dict_weights_func)
        plot_graph(X_func, size_func, file_name)
    return X_func, dict_weights_func, size_func


def pickle2pygraph(Graph, graph_type=None):
    """
    This script reads parameters about a generated pickle graph (Graph) and creates:
    - newGraph_list, -- a list containing the cc of Graph
    - mapping_list,  -- f: [1,..., len(Graph(cc))] ----> V(Graph(cc))
    - inv_mapping_list, -- f:V(Graph(cc) ----> [1,..., len(Graph(cc))]
    - [components] -- a list containing the nodes of the connected components (eg, [{1,2,3}, {4,5}])
    """

    # defining the paths where to look for the Graph
    """
    t_string = '%.0E' % decimal.Decimal(str(t))
    if os.getcwd().split('/')[-1] == 'simplifications':
        path_ = './' + folder_name[2:] + '/' + funct + '/' + funct + '_t' + t_string + '_graph' + str(
            graph_type) + '_wm' + str(weighting_method) + '.dat'
    else:
        path_ = '../muffe_sparse_optimization/simplifications/' + folder_name[
                                                                  2:] + '/' + funct + '/' + funct + '_t' + t_string + '_graph' + str(
            graph_type) + '_wm' + str(weighting_method) + '.dat'

    # the graph

    with open(path_, 'rb') as file:
        Graph = pkl.load(file)
    """
    # n_edges=len(Graph.edges()) (deleteme)

    # Defining components

    components = list(nx.connected_components(Graph))
    components.sort(key=len, reverse=True)
    print("len(comp)", len(components))

    # Creating a graph (dict) to store the different components
    Graph_ = {}

    # Create mapping/inverse mapping (dict) to store the different mappings
    mapping = {}
    inv_mapping = {}

    # Create new_Graph to store the new relabeled graphs
    newGraph = {}

    # Defining the subgraphs induced by the elements of 'components'.
    for i in range(1, len(components) + 1):
        Graph_[i] = Graph.subgraph(list(components[i - 1]))
        newGraph[i] = {}
        mapping[i] = {}
        inv_mapping[i] = {}
        N = 1
        # cycling over the nodes of the cc component as labeled in the original Graph
        for node in Graph_[i].nodes():
            # the node 'N' (\in [1,...,len(Graph_[i].nodes())]) is mapped into the 'node' in Graph_[i] thru 'mapping'
            mapping[i][N] = node
            # the other way around thru 'inv_mapping'
            inv_mapping[i][node] = N
            N += 1
        # Defining the newGraph (i.e. graph of one cc)
        newGraph[i] = nx.Graph()
        for node in Graph_[i].nodes():
            # init weights of the nodes to be 0
            weight = 0.0
            if graph_type != "3" and graph_type != "4":
                weight = Graph_[i].nodes[node]["weight"]
            pos = Graph_[i].nodes[node]["pos"]
            # adding the 'node' in Graph with a new label given by inv_mapping to newGraph. Thus, V(newGraph)C[1,...,K].
            newGraph[i].add_node(str(inv_mapping[i][node]), weight=weight, pos=pos)
        # adding the edge in Graph to newGraph with a the labels defined for the nodes
        for edge in Graph_[i].edges():
            # weights of the edges to be the same as the ones in the original graph
            weight = Graph_[i].edges[edge[0], edge[1]]["weight"]
            newGraph[i].add_edge(
                str(inv_mapping[i][edge[0]]),
                str(inv_mapping[i][edge[1]]),
                weight=weight,
            )

    # Defining lists of outputs
    newGraph_list = [newGraph[i] for i in newGraph.keys()]
    mapping_list = [mapping[i] for i in mapping.keys()]
    inv_mapping_list = [inv_mapping[i] for i in inv_mapping.keys()]

    print("number of cc:", len(newGraph_list))
    return newGraph_list, mapping_list, inv_mapping_list, [components]


def dat2pygraph(Graph, folder_name, edge_mapping, min_, BP_weights):
    """
    This script takes dat files and convert them into a python graph.
    Main usage: translating solutions of the discrete DMK solver into python graphs
    """
    # Defining the type of weights for the output python graph: opt_tdens or opt_flux

    folder_name = "./simplifications/" + folder_name[2:]
    print("we are at", folder_name, "when executing dat2pygraph")
    if BP_weights == "BPtdens":
        opt_tdens = extracting_weights(folder_name, "output/result/opt_tdens.dat")
        dict_weights, _ = weight2dict(opt_tdens, "output/result/opt_tdens.dat")
    elif BP_weights == "BPflux":
        opt_flux = extracting_weights(folder_name, "output/result/opt_flux.dat")
        dict_weights, _ = weight2dict(opt_flux, "output/result/opt_flux.dat")
    else:
        print("BP_weights not defined.")
    opt_pot = extracting_weights(folder_name, "output/result/opt_pot.dat")
    dict_weights_pot, _ = weight2dict(opt_pot, "output/result/opt_pot.dat")
    # print(dict_weights_pot)
    # print(edge_mapping)
    # Deleting edges whose weights are not greater than max*threshold

    max_ = max(dict_weights.values())
    G_simplification = nx.Graph()
    G_simplification.add_nodes_from(Graph.nodes(data=True))
    # print(G_simplification.nodes())
    for key in dict_weights.keys():
        weight = dict_weights[key]
        if weight > min_ * max_:
            G_simplification.add_edge(
                edge_mapping[key][0], edge_mapping[key][1], weight=weight
            )
            potential = dict_weights_pot[float(edge_mapping[key][0])]
            G_simplification.nodes[edge_mapping[key][0]]["op_pot"] = potential
            potential = dict_weights_pot[float(edge_mapping[key][1])]
            G_simplification.nodes[edge_mapping[key][1]]["op_pot"] = potential

    return G_simplification


print("Now here", os.getcwd())


def pygraph2dat(
    G, sources, sinks, component_index, folder_name, mapping, input_flag=None
):
    """
    This script receives a python graph and returns the some of the files required by the discrete solver:
    - graph.dat,
    - weight.dat (unitary weight for all the edges in the graph),
    - t_dens0.dat: tdens for the edges in the graph (initial conductivity),
    - rhs.dat
    Inputs:
    - Graph,
    - sources and sinks for the subGraph Graphh(cc),
    - cc index,
    - folder_name, -- (path of the graph)
    - mapping, -- f: [1,..., len(Graph(cc))] ----> V(Graph(cc))
    - inv_mapping_list, -- f:V(Graph(cc) ----> [1,..., len(Graph(cc))]
    """
    # Creating component_indx/input folder to store the inputs for the discrete DMK
    os.chdir("./simplifications/")
    print("You are here: %s" % os.getcwd())
    new_dir0 = (
        "./runs/"
        + folder_name[2:].split("/")[-2]
        + "/"
        + folder_name[2:].split("/")[-1]
        + "/component"
        + str(component_index)
    )

    try:
        os.mkdir(new_dir0)
    except OSError:
        print("Creation of the directory %s failed." % new_dir0)
    new_dirn = new_dir0 + "/input"

    try:
        os.mkdir(new_dirn)
    except OSError:
        print("Creation of the directory %s failed." % new_dirn)

    ##############################################################
    print("Creating the .dat files for this component")

    f = open(new_dirn + "/graph.dat", "w+")
    f2 = open(new_dirn + "/tdens0.dat", "w+")
    f3 = open(new_dirn + "/weight.dat", "w+")
    ##############################################################

    if input_flag == None:
        # Moving the graph_cell.dat to the folder component_index/input
        print(
            "cp  ./"
            + folder_name.split("/")[1]
            + "/"
            + folder_name.split("/")[-2]
            + "/input/graph_cell.dat"
        )
        print(new_dirn + "/" + "graph_cell.dat")
        os.system(
            "cp  ./"
            + folder_name.split("/")[1]
            + "/"
            + folder_name.split("/")[-2]
            + "/input/graph_cell.dat"
            + "  "
            + new_dirn
            + "/graph_cell.dat"
        )
        file1_ = open(new_dirn + "/graph_cell.dat", "r")
        graph_coord_triang = file1_.readlines()
        file1_.close()

        n_nodes = int(graph_coord_triang[0][:12])
        # Reading the graph_cell.dat and copying the first part (the coordinates) into graph.dat
        line_dict = {}
        number_of_nodes = len(G.nodes())
        number_of_edges = len(G.edges())
        f.write(str(number_of_nodes) + "\n")
        aux_white = " " * (12 - len(str(number_of_edges)))
        f.write(aux_white + str(number_of_edges) + "     " + graph_coord_triang[1][12:])
        counter = 1
        for line in graph_coord_triang[2 : 2 + n_nodes]:
            line_dict[str(counter)] = line
            counter += 1
        for i in G.nodes:
            f.write(line_dict[mapping[int(i)]])

        f2.write(
            "1" + (7 - len(str(number_of_edges))) * " " + str(number_of_edges) + "\n"
        )
        f2.write("time    0.0" + "\n")
        f2.write(str(number_of_edges) + "\n")

        f3.write(
            12 * " "
            + "1"
            + (12 - len(str(number_of_edges))) * " "
            + str(number_of_edges)
            + "\n"
        )
        f3.write("time 0e30" + "\n")
        f3.write((12 - len(str(number_of_edges))) * " " + str(number_of_edges) + "\n")

    elif input_flag == "image":
        line_dict = {}
        number_of_nodes = len(G.nodes())
        number_of_edges = len(G.edges())
        f.write(str(number_of_nodes) + "\n")
        aux_white = " " * (12 - len(str(number_of_edges)))
        f.write(aux_white + str(number_of_edges) + "     " + str(2) + "\n")
        counter = 1
        for node in G.nodes():
            line = G.nodes[node]["pos"]
            line_dict[counter] = line
            counter += 1
        # print(mapping[i])
        # print(G.nodes())
        for n in G.nodes:
            # print(n,mapping[i][int(n)])
            f.write(
                str(line_dict[int(n)][0]) + "   " + str(line_dict[int(n)][1]) + "\n"
            )

        f2.write(
            "1" + (7 - len(str(number_of_edges))) * " " + str(number_of_edges) + "\n"
        )
        f2.write("time    0.0" + "\n")
        f2.write(str(number_of_edges) + "\n")

        f3.write(
            12 * " "
            + "1"
            + (12 - len(str(number_of_edges))) * " "
            + str(number_of_edges)
            + "\n"
        )
        f3.write("time 0e30" + "\n")
        f3.write((12 - len(str(number_of_edges))) * " " + str(number_of_edges) + "\n")
    """
    Getting the edge relation with the new labeling,
    Defining the tdens0 and the weights based on the index of the edges whose nodes are already relabeled
    """
    edge_mapping = {}
    visited_edges = []
    white_space = "     "
    num_of_edge = 1
    for node in G.nodes():
        for neig in G.neighbors(node):
            pair1 = (neig, node)
            pair2 = (node, neig)
            if pair1 not in visited_edges or pair2 not in visited_edges:
                # Defining edge_mapping: this takes two numbers in [1,...,len(Graph(cc).edges)] and outputs the edge on Graph(cc) s.t.
                # the labelings are coherent with inv_mapping funct
                edge_mapping[num_of_edge] = [node, neig]
                f.write(
                    "  " + node + white_space + neig + "  " + str(num_of_edge) + "\n"
                )
                tdens = str(G[node][neig]["weight"])
                x = Decimal(tdens)
                y = "{:.9E}".format(x)
                f2.write(str(num_of_edge) + " " + y + "\n")
                f3.write(
                    (12 - len(str(num_of_edge))) * " "
                    + str(num_of_edge)
                    + "   "
                    + "1.00E-000"
                    + "\n"
                )
                num_of_edge += 1
                visited_edges.append((neig, node))
                visited_edges.append((node, neig))
    f2.write("time 1e+30")
    f3.write("time 1e30")
    f.close()
    f2.close()
    f3.close()

    # Creating rhs.dat: 1/(number of sources (sinks)) for each element in the source (sink) list
    str0 = 11 * " "
    f = open(new_dirn + "/rhs.dat", "w+")
    num_sources_sinks = str(len(sources) + len(sinks))
    print("G", len(G.nodes()))
    print("num ss", num_sources_sinks)
    f.write(str0 + "1" + str0[:8] + str(len(G.nodes())) + "\n")
    f.write(" time -1e30 \n")
    f.write((12 - len(num_sources_sinks)) * " " + num_sources_sinks + "\n")
    for node in sorted(list(set(sources).union(set(sinks)))):
        if node in sources:
            f.write(
                (12 - len(str(node))) * " "
                + str(node)
                + (21 - len(str(1 / len(sources)))) * " "
                + "  "
                + str(1 / len(sources))
                + "\n"
            )
        elif node in sinks:
            f.write(
                (12 - len(str(node))) * " "
                + str(node)
                + (21 - len(str(-1 / len(sinks)))) * " "
                + " -"
                + str(1 / len(sinks))
                + "\n"
            )
    f.write(" time 1e30")
    f.close()
    return edge_mapping


def fixing_weight_file(file_name):
    """
    This script adds the 'time 1e30' at the end of a dat file
    """
    # read the file into a list of lines

    lines = open(file_name, "r").readlines()

    # now edit the last line of the list of lines

    new_last_line = " time 1e30"
    lines[-1] = new_last_line

    # now write the modified list back out to the file

    open(file_name, "w").writelines(lines)


def using_graph2incidence_matrix(folder_name, index, weight_flag=None):
    """
    This script makes the files needed to use the discrete DMK solver (muffe_sparse_opt:
    graph.dat, matrix.dat, weight.dat, kernel.dat, length.dat
    """

    # main script for the generation of the files

    program = "../geometry/graph2incidence_matrix/graph2incidence_matrix.out"

    # folder path

    folder_name_here = (
        "./simplifications/runs/"
        + folder_name[2:].split("/")[-2]
        + "/"
        + folder_name[2:].split("/")[-1]
    )

    # defining the paths of the files

    graph = folder_name_here + "/component" + str(index) + "/input" + "/graph.dat"
    matrix = folder_name_here + "/component" + str(index) + "/input" + "/matrix.dat"

    # (deciding whether the length.dat will be called weight.dat or length.dat)

    if weight_flag == "length":
        length = folder_name_here + "/component" + str(index) + "/input" + "/weight.dat"
    else:
        length = folder_name_here + "/component" + str(index) + "/input" + "/length.dat"

    kernel = folder_name_here + "/component" + str(index) + "/input" + "/kernel.dat"

    # generating the files defined before

    command = (
        str(program)
        + " "
        + str(graph)
        + " "
        + str(matrix)
        + " "
        + str(length)
        + " "
        + str(kernel)
    )
    os.system(command)

    # fixing the file

    if weight_flag == "length":
        fixing_weight_file(length)


def updating_beta_discrete(beta):
    """
    This script generates a new pflux.dat file in the cwd = Tests/
    """
    print("defining pflux.dat for beta=", beta)

    # Writing the file

    f = open("./simplifications/par_files/pflux.dat", "w+")
    f.write("1  1 \n")
    f.write("time  0.0 \n")
    f.write("1 !ninputs \n")
    beta_str = str(beta)
    f.write("1 " + beta_str + "\n")
    f.write("time 1.e30 \n")
