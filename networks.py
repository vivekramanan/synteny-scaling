"""
Network Plotting Functions
"""
import pandas as pd
import networkx as nx
import numpy as np
import plotly.graph_objects as go
import plotly.express as px

def plotNetwork(G, k, taxColor, taxD, labels, title):
    edge_x = []
    edge_y = []
    for edge in G.edges():
        x0, y0 = k[edge[0]]
        x1, y1 = k[edge[1]]
        edge_x.append(x0)
        edge_x.append(x1)
        edge_x.append(None)
        edge_y.append(y0)
        edge_y.append(y1)
        edge_y.append(None)

    edge_trace = go.Scatter(
        x=edge_x, y=edge_y,
        line=dict(width=0.5, color='#888'),
        hoverinfo='none',
        mode='lines')

    node_x = []
    node_y = []
    for node in G.nodes():
        x, y = k[node]
        node_x.append(x)
        node_y.append(y)

    node_trace = go.Scatter(
        x=node_x, y=node_y,
        mode='markers',
        hoverinfo='text',
        marker=dict(
            size=10,
            line_width=2,
            color=list(taxColor.values())))

    node_tax= []
    node_adjacencies = []
    node_text = []
    for node, adjacencies in enumerate(G.adjacency()):
        node_tax.append(taxColor[taxD[labels[node]]])
        node_adjacencies.append(len(adjacencies[1]))
        node_text.append(str(labels[node])+': '+str(len(adjacencies[1])))

    node_trace.marker.color = node_tax
    #node_trace.marker.size = node_adjacencies
    node_trace.text = node_text

    fig = go.Figure(data=[edge_trace, node_trace],
                 layout=go.Layout(
                    title='<br>'+title,
                    titlefont_size=16,
                    showlegend=False,
                    hovermode='closest',
                    margin=dict(b=20,l=5,r=5,t=50),
                    paper_bgcolor='rgba(0,0,0,0)',
                    plot_bgcolor='rgba(0,0,0,0)',
                    annotations=[ dict(
                        text="k=10, Kamada-Kawai Layout",
                        showarrow=False,
                        xref="paper", yref="paper",
                        x=0.005, y=-0.002 ) ],
                    xaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
                    yaxis=dict(showgrid=False, zeroline=False, showticklabels=False))
                    )
    fig.show()

def plotIsoNetwork(G, k, taxColor, taxD, labels, title):
    edge_x = []
    edge_y = []
    for edge in G.edges():
        x0, y0 = k[edge[0]]
        x1, y1 = k[edge[1]]
        edge_x.append(x0)
        edge_x.append(x1)
        edge_x.append(None)
        edge_y.append(y0)
        edge_y.append(y1)
        edge_y.append(None)

    edge_trace = go.Scatter(
        x=edge_x, y=edge_y,
        line=dict(width=0.5, color='#888'),
        hoverinfo='none',
        mode='lines')

    node_x = []
    node_y = []
    for node in G.nodes():
        x, y = k[node]
        node_x.append(x)
        node_y.append(y)

    node_trace = go.Scatter(
        x=node_x, y=node_y,
        mode='markers',
        hoverinfo='text',
        marker=dict(
            size=10,
            line_width=2,
            color=list(taxColor.values())))

    node_tax= []
    node_adjacencies = []
    node_text = []
    for node, adjacencies in enumerate(G.adjacency()):
        node_tax.append(taxColor[taxD[labels[adjacencies[0]]]])
        node_adjacencies.append(len(adjacencies[1]))
        node_text.append(str(labels[adjacencies[0]])+': '+str(len(adjacencies[1])))

    node_trace.marker.color = node_tax
    #node_trace.marker.size = node_adjacencies
    node_trace.text = node_text

    fig = go.Figure(data=[edge_trace, node_trace],
                 layout=go.Layout(
                    title='<br>'+title,
                    titlefont_size=16,
                    showlegend=False,
                    hovermode='closest',
                    margin=dict(b=20,l=5,r=5,t=50),
                    paper_bgcolor='rgba(0,0,0,0)',
                    plot_bgcolor='rgba(0,0,0,0)',
                    annotations=[ dict(
                        text="Spring Layout",
                        showarrow=False,
                        xref="paper", yref="paper",
                        x=0.005, y=-0.002 ) ],
                    xaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
                    yaxis=dict(showgrid=False, zeroline=False, showticklabels=False))
                    )
    fig.show()