#!/usr/bin/env python

import plotly.graph_objs as go

"""

# Add shapes
fig2.layout.shapes=[
        go.layout.Shape(
            type="rect",
            x0=1/8,
            y0=1/4,
            x1=3/8,
            y1=3/4,
            line=dict(
                color="RoyalBlue",
            ),
        ),
        # filled Rectangle
        go.layout.Shape(
            type="rect",
            x0=5/8,
            y0=1/4,
            x1=7/8,
            y1=3/4,
            line=dict(
                color="RoyalBlue",
            ),
        ),
    ]

iplot(fig2, filename='networkx')
#print("The number of out nodes is: ", int(out))
#fig.show()  

"""

"""
X_dir = directing(newGraph)

source_nodes_2=[]
sink_nodes_2=[]
nodes2=list(X_dir.nodes())

"""
##########################################################################


# if __name__ == "__main__" :
x1, y1, x2, y2 = 1.0 / 8.0, 1 / 4, 3 / 8, 3 / 4
x3, y3, x4, y4 = 5.0 / 8.0, 1 / 4, 7 / 8, 3 / 4


def Source(x1, y1, x2, y2, x, y):
    if x > x1 and x < x2 and y > y1 and y < y2:
        return True
    else:
        return False


def Sink(x3, y3, x4, y4, x, y):
    if x > x3 and x < x4 and y > y3 and y < y4:
        return True
    else:
        return False


def getting_sources_sinks(Graph):
    print("getting_sources_sinks. Version: 1.2")
    source_nodes = []
    sink_nodes = []
    for node in Graph.nodes():
        x, y = Graph.node[node]["pos"]
        # node_trace_1['x'] += tuple([x])
        # node_trace_1['y'] += tuple([y])

        if (Graph.out_degree(node) > 0 and Graph.in_degree(node) == 0) or (
            Graph.out_degree(node) == 0 and Graph.in_degree(node) > 0
        ):
            # print(x, y, "node in the source")
            if Source(x1, y1, x2, y2, x, y):
                source_nodes.append(node)
                # print(node, "Node in the source, and this is a beggining")
            elif Sink(x3, y3, x4, y4, x, y):
                # print(x, y, "Node in the sink!")
                sink_nodes.append(node)
                # print(node, "Node in the sink, and this is the end")
    # Add shapes
    source_trace = go.Scatter(
        x=[], y=[], name="source", line=dict(color="RoyalBlue"), mode="lines"
    )

    for y_ in [1 / 4]:
        for x_ in [1 / 8, 3 / 8]:
            source_trace["x"] += tuple([x_])
            source_trace["y"] += tuple([y_])
    for y_ in [3 / 4]:
        for x_ in [3 / 8, 1 / 8]:
            source_trace["x"] += tuple([x_])
            source_trace["y"] += tuple([y_])
    source_trace["x"] += tuple([1 / 8])
    source_trace["y"] += tuple([1 / 4])
    # filled Rectangle
    sink_trace = go.Scatter(
        x=[], y=[], name="sink", line=dict(color="RoyalBlue"), mode="lines"
    )

    for y_ in [1 / 4]:
        for x_ in [5 / 8, 7 / 8]:
            sink_trace["x"] += tuple([x_])
            sink_trace["y"] += tuple([y_])

    for y_ in [3 / 4]:
        for x_ in [7 / 8, 5 / 8]:
            sink_trace["x"] += tuple([x_])
            sink_trace["y"] += tuple([y_])
    sink_trace["x"] += tuple([5 / 8])
    sink_trace["y"] += tuple([1 / 4])
    return source_nodes, sink_nodes, source_trace, sink_trace
