import sys
from itertools import takewhile
from perseuspy import pd, nx, write_networks
import numpy as np
import perseuspy.parameters as params
import phos.data.network as networker
import phos.algo.anat as anat

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='PHOTON: PHOsphoproteomic dissecTiOn using Networks')
    parser.add_argument('paramfile', type=argparse.FileType('r'))
    parser.add_argument('infile', type=argparse.FileType('r'))
    parser.add_argument('outfolder', type=str)
    args = parser.parse_args()
    print('reading parameters')
    paramFile = params.parse_parameters(args.paramfile)
    parameters = {
            'ppi-network' : {
                'ppi_network' : params.fileParam(paramFile, "Network"),
                'confidence' : params.doubleParam(paramFile, 'Network confidence'),
                'degree_threshold' : params.intParam(paramFile, 'Network node degree threshold')
                }}
    data = pd.read_perseus(args.infile)
    geneid_column = params.singleChoiceParam(paramFile, 'GeneID')
    if type(geneid_column) is int and geneid_column < 0:
        print("GeneID column was not chosen")
        sys.exit(1)
    data = data[pd.notnull(data[geneid_column])]
    try:
        data[geneid_column] = data[geneid_column].astype(int)
    except ValueError:
        print('\n'.join(["GeneID column could not be converted to integer values",
            "Please remove empty, or non-numeric values such as ';'"]))
        sys.exit(1)
    anchor = params.stringParam(paramFile, 'Signaling source')
    if anchor is not None:
        try:
            anchor = int(anchor)
        except ValueError:
            print("Cannot convert signaling source {} to entrez gene id. Please enter a number, or leave blank".format(anchor))
            sys.exit(1)
    sig_column, subparams = params.singleChoiceWithSubParams(paramFile, 'Select terminals from')
    if type(sig_column) is int and sig_column < 0:
        print("Please select a categorical column for terminals selection.")
        sys.exit(1)
    sig_values = params.multiChoiceParam(subparams, 'Select')
    if len(sig_values) < 1:
        print("Please select at least one value for choosing terminals.")
        sys.exit(1)
    terminals = data[data[sig_column].isin(sig_values)][geneid_column]
    if len(terminals) < 1:
        print("No terminals were found. Please make sure that the selected values are present in the table.")
        sys.exit(1)
    network_undirected = networker.load(**parameters['ppi-network'])
    subnet = anat.remote_network('Perseus', network_undirected, terminals, anchor = anchor)
    subnet.to_csv('subnet.csv')
    G = nx.from_pandas_dataframe(subnet, 's', 't', edge_attr=True)
    for node in G:
        if node in terminals:
            node_type = 'terminal'
        elif node == anchor:
            node_type = 'anchor'
        else:
            node_type = 'connector'
        G.node[node]['Type'] = node_type
    networks_table, networks = nx.to_perseus([G])
    print('writing results')
    write_networks(args.outfolder, networks_table, networks)
