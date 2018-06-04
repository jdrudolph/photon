import sys
from itertools import takewhile
from perseuspy import pd
import numpy as np
import perseuspy.parameters as params

try:
    import phos.data.network as networker
    import phos.algo.activity as activity
except ImportError:
    print("Could not import the photon python package")
    sys.exit(1)

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='PHOTON: PHOsphoproteomic dissecTiOn using Networks')
    parser.add_argument('paramfile', type=argparse.FileType('r'))
    parser.add_argument('infile', type=argparse.FileType('r'))
    parser.add_argument('outfile', type=argparse.FileType('w'))
    args = parser.parse_args()
    print('reading parameters')
    paramFile = params.parse_parameters(args.paramfile)
    parameters = {
            'ppi-network' : {
                'ppi_network' : params.fileParam(paramFile, "Network"),
                'confidence' : params.doubleParam(paramFile, 'Network confidence'),
                'degree_threshold' : params.intParam(paramFile, 'Network node degree threshold')
                },
            'activity' : {
                'min_size' : params.intParam(paramFile, 'Required number of observations'),
                'permutations' : params.intParam(paramFile, 'Number of permutations')
                }}
    data = pd.read_perseus(args.infile)
    column_names = data.columns.get_level_values('Column Name') 
    data_columns = [column_names[i] for i, t in takewhile(lambda tup: tup[1] is np.dtype('float64'), enumerate(data.dtypes))]
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
    aa_column = params.singleChoiceParam(paramFile, 'Amino acid')
    if type(aa_column) is int and aa_column < 0:
        print("Amino acid column was not chosen")
        sys.exit(1)
    pos_column = params.singleChoiceParam(paramFile, 'Position')
    if  type(pos_column) is int and pos_column < 0:
        print("Position column was not chosen")
        sys.exit(1)
    rename = {geneid_column : 'GeneID', aa_column : 'Amino.Acid', pos_column : 'Position'}
    columns = ['GeneID', 'Amino.Acid', 'Position', 'avg']
    results = []
    for data_column in data_columns:
        print('calculating scores for', data_column)
        rename[data_column] = 'avg'
        exp = data[[geneid_column, aa_column, pos_column, data_column]].rename(columns = rename).dropna()
        exp.columns = exp.columns.get_level_values('Column Name')
        network_undirected = networker.load(**parameters['ppi-network'])
        network = networker.to_directed(network_undirected)
        scores = activity.empiric(exp, network, **parameters['activity'])
        scores['Column Name'] = data_column
        results.append(scores)
    _result = pd.concat(results)
    result = _result.pivot_table('score_empiric', 'GeneID', 'Column Name')
    for num_column in ['t', 'p_greater', 'p_lesser', 'p_twosided', 'q_greater', 'q_lesser', 'q_twosided']:
        values = (_result
                .pivot_table(num_column, 'GeneID', 'Column Name')
                .rename(columns = lambda x : '{} {}'.format(x, num_column))
                .apply(lambda col: col.astype(float)))
        result = result.join(values)
    for cat_column in ['rej_greater', 'rej_lesser', 'rej_twosided', 'Significant']:
        significant = (_result
                .pivot_table(cat_column, 'GeneID', 'Column Name')
                .rename(columns = lambda x : '{} {}'.format(x, cat_column))
                .fillna('Insufficient data')
                .apply(lambda col: col.astype('category')))
        result = result.join(significant)
    result = result.reset_index();
    result['GeneID'] = result['GeneID'].astype(int).astype(str)
    print('writing results')
    result.to_perseus(args.outfile, main_columns=data_columns)
