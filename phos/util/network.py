import pandas as pd

def add_edge_is_directed(network):
    """
    Determine the edge directions in an undirected network
    by adding an 'is_directed' column to the network.
    """
    reverse_network = network[['kin', 'sub']].rename(columns= {'kin': 'sub', 'sub': 'kin'})
    network = network.merge(reverse_network, on=['kin', 'sub'], how='left', indicator=True)
    directed = network[network['_merge'] == 'left_only'].drop('_merge', 1)
    undirected = network[network['_merge'] == 'both'].drop('_merge', 1)
    undirected_collapsed = (undirected[undirected['kin'] <= undirected['sub']]
            .drop_duplicates(subset=['kin', 'sub']))
    undirected_collapsed['is_directed'] = False
    directed['is_directed'] = True
    return pd.concat([undirected_collapsed, directed])

