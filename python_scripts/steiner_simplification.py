import networkx as nx
from networkx.algorithms.approximation import steiner_tree
import pickle as pkl
import numpy as np
import matplotlib.pyplot as plt
import plotly.tools as tls

# import chart_studio.plotly as py
# import plotly.express as px
import plotly.graph_objs as go

# from graph2plotly import *
import copy
from Getting_sources_and_sinks import *


def directing(G):
    """
    Script to make directed copies of the optimal graphs. 
    The direction is taken to be consistent with the from-source-to-sink dynamics.
    """
    G_directed = nx.DiGraph()
    G_directed.add_nodes_from(G.nodes(data=True))
    for edge in G.edges():
        x0, y0 = G.node[edge[0]]["pos"]
        x1, y1 = G.node[edge[1]]["pos"]
        if x0 < x1:
            G_directed.add_edge(edge[0], edge[1])
        else:
            G_directed.add_edge(edge[1], edge[0])
    return G_directed


def fixing_weights(G):
    max_ = max([G.edges[edge]["weight"] for edge in G.edges()])
    print("Maximum weight:", max_)
    for edge in G.edges():
        G.edges[edge]["weight"] = 1 - G.edges[edge]["weight"] / max_
        # print(G.edges[edge]['weight'])
    return G


# def graph_dep_threshold(G):


def steiner_simpl(G, plotting):
    """
    Algorithm used to simplify the optimal graphs via Steiner tree technique
    (minimum spanning tree and shortest path)
    """
    G = fixing_weights(G)
    terminal_selector = "betweenness_centrality"  # 'directed_graph', 'betweenness_centrality', 'degree_centrality''
    print(
        "Chosen criteria for the classification of the sources and sinks:",
        terminal_selector,
    )
    print("directing G")
    components = list(nx.connected_components(G))
    components.sort(key=len, reverse=True)
    print("size of the graph:", len(G.nodes()))
    print("number of connected components:", len(components))
    node_trace_source_sink = go.Scatter(
        x=[],
        y=[],
        mode="markers",
        name="source/sink",
        text="",
        hoverinfo="text",
        marker=dict(color="green", size=10, line=dict(width=2)),
    )
    for i in range(len(components)):
        print(i, "-component has length equal to:", len(components[i]))
    if len(components) == 1:
        component1 = components[0]
        if terminal_selector == "directed_graph":
            G_directed = directing(G)
            G_directed_1 = G_directed.copy()
            (
                source_nodes_1,
                sink_nodes_1,
                source_trace,
                sink_trace,
            ) = getting_sources_sinks(G_directed_1)
            print("source nodes:", source_nodes_1)
            print("sink nodes:", sink_nodes_1)
            terminal_nodes_1 = source_nodes_1 + sink_nodes_1
            ss_nodes = source_nodes_1 + sink_nodes_1
            print("Steiner tree 1 cc")
            G_simplification = steiner_tree(
                G_directed_1.to_undirected(), terminal_nodes_1
            )
        elif terminal_selector != "directed_graph":
            G_1 = G.copy()
            if terminal_selector == "betweenness_centrality":
                bn1 = nx.betweenness_centrality(G_1, normalized=True)
                min_bn1 = min(bn1.values())
                kind_of_leaf_nodes_1 = [
                    key for key in bn1.keys() if bn1[key] == min_bn1
                ]
                sub1 = G_1.subgraph(kind_of_leaf_nodes_1)
                edges_to_remove_1 = list(sub1.edges())
                print("edges in the graph of sources 1:", sub1.edges())
                G_1.remove_edges_from(edges_to_remove_1)
                G.remove_edges_from(edges_to_remove_1)
            elif terminal_selector == "degree_centrality":
                dg1 = nx.degree_centrality(G_1)
                min_dg1 = min(dg1.values())
                kind_of_leaf_nodes_1 = [
                    key for key in dg1.keys() if dg1[key] == min_dg1
                ]

            source_nodes_1 = []
            sink_nodes_1 = []
            for key in kind_of_leaf_nodes_1:
                x, y = G.node[key]["pos"]
                if Source(x1, y1, x2, y2, x, y):
                    source_nodes_1.append(key)
                    # print(node, "Node in the source, and this is a beggining")
                elif Sink(x3, y3, x4, y4, x, y):
                    # print(x, y, "Node in the sink!")
                    sink_nodes_1.append(key)
            terminal_nodes_1 = source_nodes_1 + sink_nodes_1

            G_simplification = steiner_tree(G_1, terminal_nodes_1)

            ss_nodes = terminal_nodes_1
        _, _, source_trace, sink_trace = getting_sources_sinks(nx.DiGraph())
    elif len(components) > 1:
        component1 = components[0]
        component2 = components[1]
        rest_of_components = []
        for comp in components[2:]:
            rest_of_components = rest_of_components + list(comp)
        if terminal_selector == "directed_graph":
            G_directed = directing(G)
            G_directed_1 = G_directed.copy()
            G_directed_1.remove_nodes_from(list(component2) + rest_of_components)
            # Lower component
            G_directed_2 = G_directed.copy()
            G_directed_2.remove_nodes_from(list(component1) + rest_of_components)

            (
                source_nodes_1,
                sink_nodes_1,
                source_trace,
                sink_trace,
            ) = getting_sources_sinks(G_directed_1)
            source_nodes_2, sink_nodes_2, _, _ = getting_sources_sinks(G_directed_2)
            print("source nodes C1", source_nodes_1)
            print("sink nodes C1", sink_nodes_1)
            print("source nodes C2", source_nodes_2)
            print("sink nodes C2", sink_nodes_2)
            terminal_nodes_1 = source_nodes_1 + sink_nodes_1
            terminal_nodes_2 = source_nodes_2 + sink_nodes_2
            ss_nodes = source_nodes_1 + source_nodes_2 + sink_nodes_1 + sink_nodes_2

            # Steiner algorithm
            print("steiner tree algorithm in the first component")
            if source_nodes_1 != [] and sink_nodes_1 != []:
                G_1_simplification = steiner_tree(
                    G_directed_1.to_undirected(), terminal_nodes_1
                )
            else:
                print(
                    "this connected component does NOT connect the source and the sink. Skipping steiner simplification"
                )
                G_1_simplification = nx.Graph()
            # print(G_1_simplification.nodes(data=True))
            print("steiner tree algorithm in the second component")
            if source_nodes_2 != [] and sink_nodes_2 != []:
                G_2_simplification = steiner_tree(
                    G_directed_2.to_undirected(), terminal_nodes_2
                )
            else:
                print(
                    "this connected component does NOT connect the source and the sink. Skipping steiner simplification"
                )
                G_2_simplification = nx.Graph()
        elif terminal_selector != "directed_graph":
            G_1 = G.copy()
            G_1.remove_nodes_from(list(component2) + rest_of_components)
            G_2 = G.copy()
            G_2.remove_nodes_from(list(component1) + rest_of_components)
            if terminal_selector == "betweenness_centrality":
                N = 1
                print('threshold for "kind-of-leaf" nodes:', N)
                bn1 = nx.betweenness_centrality(G_1, normalized=True)
                bn2 = nx.betweenness_centrality(G_2, normalized=True)
                bn1_list = list(set(bn1.values()))
                bn1_list.sort()
                print(bn1_list[0], bn1_list[1])
                bn2_list = list(set(bn2.values()))
                bn2_list.sort()

                kind_of_leaf_nodes_1 = [
                    key for key in bn1.keys() if bn1[key] in bn1_list[0:N]
                ]  # or bn1[key]==bn1_list[1]) ]# or bn1[key]==bn1_list[2])]
                kind_of_leaf_nodes_2 = [
                    key for key in bn2.keys() if bn2[key] in bn2_list[0:N]
                ]  # or bn2[key]==bn2_list[1]) ]#or bn2[key]==bn2_list[2])]
                """
                kind_of_leaf_nodes_1=[]
                for key in bn1.keys():
                    print(min_bn1,1.1*min_bn1)

                    if (bn1[key]>=min_bn1 and bn1[key]<1.1*min_bn1):
                        kind_of_leaf_nodes_1.append(key)
                kind_of_leaf_nodes_2=[]
                for key in bn2.keys():
                    if (bn2[key]>=min_bn2 and bn2[key]<1.1*min_bn2):
                        kind_of_leaf_nodes_2.append(key)
                """
                # ...=[key for key in bn1.keys() if bn1[key]>=min_bn1 and bn1[key]<1.1*min_bn1]
                # ...=[key for key in bn2.keys() if (bn2[key]>=min_bn2 and bn2[key]<1.1*min_bn2)]
                print(kind_of_leaf_nodes_1)
                print(kind_of_leaf_nodes_2)
                # removing the connections between the sources and sinks
                sub1 = G_1.subgraph(kind_of_leaf_nodes_1)
                sub2 = G_2.subgraph(kind_of_leaf_nodes_2)
                edges_to_remove_1 = list(sub1.edges())
                edges_to_remove_2 = list(sub2.edges())
                print("edges in the graph of sources 1:", sub1.edges())
                print("edges in the graph of sources 2:", sub2.edges())
                G_1.remove_edges_from(edges_to_remove_1)
                G_2.remove_edges_from(edges_to_remove_2)
                G.remove_edges_from(edges_to_remove_1 + edges_to_remove_2)
            elif terminal_selector == "degree_centrality":
                dg1 = nx.degree_centrality(G_1)
                dg2 = nx.degree_centrality(G_2)
                min_dg1 = min(dg1.values())
                min_dg2 = min(dg2.values())
                kind_of_leaf_nodes_1 = [
                    key for key in dg1.keys() if dg1[key] == min_dg1
                ]
                kind_of_leaf_nodes_2 = [
                    key for key in dg2.keys() if dg2[key] == min_dg2
                ]
        source_nodes_1 = []
        source_nodes_2 = []
        sink_nodes_1 = []
        sink_nodes_2 = []
        if len(components) > 1:
            for key in kind_of_leaf_nodes_1:
                x, y = G.node[key]["pos"]
                if Source(x1, y1, x2, y2, x, y):
                    source_nodes_1.append(key)
                    # print(node, "Node in the source, and this is a beggining")
                elif Sink(x3, y3, x4, y4, x, y):
                    # print(x, y, "Node in the sink!")
                    sink_nodes_1.append(key)
            for key in kind_of_leaf_nodes_2:
                flag = False
                x, y = G.node[key]["pos"]
                if Source(x1, y1, x2, y2, x, y):
                    flag = True
                    source_nodes_2.append(key)
                    # print(node, "Node in the source, and this is a beggining")
                elif Sink(x3, y3, x4, y4, x, y):
                    flag = True
                    sink_nodes_2.append(key)
                if flag == True:
                    node_trace_source_sink["x"] += tuple([x])
                    node_trace_source_sink["y"] += tuple([y])
        terminal_nodes_1 = source_nodes_1 + sink_nodes_1
        terminal_nodes_2 = source_nodes_2 + sink_nodes_2
        G_1_simplification = steiner_tree(G_1, terminal_nodes_1)
        G_2_simplification = steiner_tree(G_2, terminal_nodes_2)
        ss_nodes = terminal_nodes_1 + terminal_nodes_2
        _, _, source_trace, sink_trace = getting_sources_sinks(nx.DiGraph())
        G_simplification = nx.union(G_1_simplification, G_2_simplification)
    for key in kind_of_leaf_nodes_1:
        flag = False
        x, y = G.node[key]["pos"]
        if Source(x1, y1, x2, y2, x, y):
            flag = True
            source_nodes_1.append(key)
        elif Sink(x3, y3, x4, y4, x, y):
            flag = True
        sink_nodes_1.append(key)
        if flag == True:
            node_trace_source_sink["x"] += tuple([x])
            node_trace_source_sink["y"] += tuple([y])

    print("saving plots")
    graph2plotly(
        folder_name,
        t,
        graph_type,
        [G, G_simplification],
        "both",
        ["blue", "red"],
        "steiner_simpl_",
        [node_trace_source_sink, source_trace, sink_trace],
    )

    return G_simplification


print("steiner_simplification. Version: 1.2")


"""
Test
"""
test = "no"
for i in range(1, 2):
    case = str(11)
    if test == "yes":
        if case == "1":
            folder_name = "../../../data/1_refinement/run_sum_parabolas_3"
        elif case == "2":
            folder_name = "../../../data/1_refinement/run_parabola_4"
        elif case == "3":
            folder_name = "../../../data/1_refinement/run_multiple_deltas_2"
        elif case == "4":
            folder_name = "../../../data/1_refinement/run_parabola_8"
        elif case == "5":
            folder_name = "../../../data/1_refinement/run_circle_2"
        elif case == "6":
            folder_name = "../../../data/1_refinement/run_parabola_9"
        elif case == "7":
            folder_name = "../../../data/1_refinement/run_parabola_2"
        elif case == "8":
            folder_name = "../../../data/1_refinement/run_uniform_3"
        elif case == "9":
            folder_name = "../../../data/1_refinement/run_sine"
        elif case == "10":
            folder_name = "../../../data/1_refinement/run_cosine"
        elif case == "11":
            folder_name = "../../../data/1_refinement/parabola_1_1ref"

        t = 0.01
        graph_type = "1"

        path_ = (
            folder_name
            + "/opt_tdens/"
            + "graph_pc"
            + str(int(t * 100))
            + "_g"
            + graph_type
            + ".dat"
        )
        print(path_)
        with open(path_, "rb") as file:
            Graph = pkl.load(file)
        # print(len(Graph.edges()))
        """
        components=list(nx.connected_components(Graph))
        component1 = components[0] 
        component2 = components[1]
        component3 = components[2]
        G_test=Graph.copy()
        G_test.remove_nodes_from(list(component2)+list(component3))

        print('computing steiner simplification')
        steiner_simpl(G_test,"yes")
        """
        print("computing steiner simplification")
        steiner_simpl(Graph, "yes")

testing_fixing_weights = "no"
if testing_fixing_weights == "yes":
    folder_name = "../../../data/1_refinement/parabola_1_1ref"
    t = 0.01
    graph_type = "1"
    path_ = (
        folder_name
        + "/opt_tdens/"
        + "graph_pc"
        + str(int(t * 100))
        + "_g"
        + graph_type
        + ".dat"
    )
    print(path_)
    with open(path_, "rb") as file:
        Graph = pkl.load(file)

    fixing_weights(Graph)
