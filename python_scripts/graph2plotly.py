import plotly.tools as tls

# import chart_studio.plotly as py
import plotly.graph_objs as go


def graph2plotly(
    folder_name,
    t,
    graph_type,
    listOfGraphs,
    style_,
    color_list,
    funct,
    graph_name=None,
    extra_data=None,
):
    data = []
    data_ = {}
    node_trace = {}
    edge_trace = {}
    for i in range(len(listOfGraphs)):
        if i == 0:
            name_ = "original_graph"
        elif i == 1 and graph_name[:5] == "stein":

            name_ = "steiner_simplification"
        else:
            name_ = "graph_N" + str(i)
        col = color_list[i]
        X = listOfGraphs[i].copy()

        node_trace[i] = go.Scatter(
            x=[],
            y=[],
            name=name_,
            mode="markers",
            # hoverinfo='text',
            marker=dict(color=col, size=12, line=dict(width=2)),
        )

        for node in X.nodes():
            x, y = X.node[node]["pos"]
            node_trace[i]["x"] += tuple([x])
            node_trace[i]["y"] += tuple([y])

        edge_trace[i] = go.Scatter(
            x=[],
            y=[],
            name=name_,
            line=dict(width=2, color=col),
            hoverinfo="none",
            mode="lines",
        )

        for edge in X.edges():
            x0, y0 = X.node[edge[0]]["pos"]
            x1, y1 = X.node[edge[1]]["pos"]
            edge_trace[i]["x"] += tuple([x0, x1, None])
            edge_trace[i]["y"] += tuple([y0, y1, None])
        data_[i] = [node_trace[i], edge_trace[i]]
        data.append(node_trace[i])
        data.append(edge_trace[i])

        if style_ == "separated" or style_ == "both":
            if extra_data != None:  # and i==len(listOfGraphs)-1:
                data_[i] = data_[i] + extra_data

            fig = go.Figure(
                data=data_[i],
                layout=go.Layout(
                    title=graph_name + "together_threshold " + str(t * 100) + "%",
                    titlefont=dict(size=16),
                    showlegend=True,
                    hovermode="closest",
                    autosize=False,
                    width=900,
                    height=700,
                    margin=dict(b=20, l=20, r=20, t=40),
                    xaxis=dict(showgrid=True, zeroline=True, showticklabels=True),
                    yaxis=dict(showgrid=True, zeroline=True, showticklabels=True),
                ),
            )

            fig.write_image(
                folder_name
                + "/"
                + funct
                + "/"
                + name_
                + "_pc"
                + str(int(t * 100))
                + "_graph"
                + graph_type
                + ".png"
            )  ###

    if style_ == "together" or style_ == "both":
        if extra_data != None and i == len(listOfGraphs) - 1:
            data = data + extra_data
        fig = go.Figure(
            data=data,
            layout=go.Layout(
                title=" ",  # graph_name+ 'together_threshold ' +str(t*100) +'%',
                titlefont=dict(size=16),
                showlegend=True,
                hovermode="closest",
                autosize=False,
                width=900,
                height=700,
                margin=dict(b=20, l=20, r=20, t=40),
                xaxis=dict(showgrid=True, zeroline=True, showticklabels=True),
                yaxis=dict(showgrid=True, zeroline=True, showticklabels=True),
            ),
        )
        if graph_name == None:
            aux_name = "Graph"
        else:
            aux_name = graph_name
        print("saving plot at", folder_name)
        # fig.update_layout(plot_bgcolor='rgb(255,255,255)')
        fig.write_image(
            folder_name
            + "/"
            + funct
            + "/"
            + aux_name
            + "pc"
            + str(int(t * 100))
            + "_graph"
            + graph_type
            + ".png"
        )

        # iplot(fig2, filename='networkx')


print("graph2plotly. Version: 1.2")
