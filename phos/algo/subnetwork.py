import os
import shutil
import json
import goenrich
import pandas as pd
import numpy as np

import phos.algo.anat

def inference(exp, scores, network_undirected, anat, go, task_id, db, **kwargs):
    active = scores[scores['Significant']]['GeneID']
    # LOCAL
    # subnetwork = phos.algo.anat.local_network(network_undirected,
    #    active, exp, task_id=task_id, **anat)
    subnetwork = phos.algo.anat.remote_network(task_id, network_undirected,
            active, anat.get('anchor', None)).astype(int)
    if subnetwork is None:
        print("Subnetwork did not contain any edges")
        subnetwork = pd.DataFrame({'s':[-1], 't':['-1']})
    # GO enrichment
    O = goenrich.obo.ontology(db['go']['ontology'])
    annotations = goenrich.read.gene2go(db['go']['annotations'])
    nodes = pd.DataFrame({'GeneID' :
        list(set(network_undirected['kin']) | set(network_undirected['sub']))})
    background = goenrich.tools.generate_background(annotations, nodes, 'GO_ID', 'GeneID')
    query = (set(subnetwork['s']) | set(subnetwork['t'])) - set([anat.get('anchor', None)])
    goenrich.enrich.propagate(O, background, 'scores')
    go_scores = goenrich.enrich.analyze(O, query, 'scores', **go)
    return subnetwork, go_scores

import networkx as nx
from networkx.readwrite import json_graph

def create_graph(exp, scores, network, task_id, db, anchor=None):
    """ generate result graph
    :param exp: pandas.DataFrame with columns: Symbol, GeneID, Amino.Acid, Position, avg
    :param scores: protein activity scores 
    :param network: pandas.DataFrame with columns 's' -> 't'
    :param anchor: anchor id
    """
    terminals = set(scores[scores['Significant']]['GeneID'])
    G = nx.from_edgelist([[int(x) for x in y] for y in network[['s','t']].values])
    df = (exp[['Symbol', 'GeneID', 'Amino.Acid', 'Position', 'avg']]
            [exp['GeneID'].isin(network['s']) | exp['GeneID'].isin(network['t'])])
    node_attributes = (df.groupby('GeneID')
            .apply(lambda x : [{'AA': a, 'POS': int(p), 'AVG': v, 'NUM': len(x)} for a,p,v in
                zip(x['Amino.Acid'].values, x['Position'].values, x['avg'].values)]))
    geneinfo_table = pd.read_table(db['geneinfo'], comment='#', usecols=[1, 2], header=None)
    gene_name = dict(geneinfo_table.dropna().values)
    for n in G:
        G.node[n]['residues'] = node_attributes.get(n, [{'NUM' : 0}])
        G.node[n]['type'] = 'terminal' if n in terminals else 'connector'
        G.node[n]['name'] = gene_name.get(n, n)
    if anchor in G and anchor is not None:
        G.node[anchor]['type'] = 'anchor'
    return G

def draw(exp, scores, network, task_id, template_dir, db, anchor=None):
    """ generate interactive result graph based on d3js
    
    check `phos/algo/anat.html` for the javascript part of the visualization

    :param exp: pandas.DataFrame with columns: Symbol, GeneID, Amino.Acid, Position, avg
    :param scores: protein activity scores 
    :param network: pandas.DataFrame with columns 's' -> 't'
    :param anchor: anchor id
    :returns HTML: the generated graph
    """
    G = create_graph(exp, scores, network, task_id, db, anchor)
    graph = json_graph.node_link_data(G)
    node_index = {n['id']: i for i,n in enumerate(graph['nodes'])}
    graph['links'] = [{'source': node_index[link['source']], 'target': node_index[link['target']]} for link in graph['links']]
    from jinja2 import Environment, FileSystemLoader
    env = Environment(loader=FileSystemLoader(template_dir))
    HTML = env.get_template('result.html').render(task_id=task_id, graph=json.dumps(graph))
    return HTML

def _rnd_go(exp, _score, network_undirected, anat, go, seed):
    np.random.seed(seed)
    rnd = (_score[['Significant']]
            .reindex(np.random.permutation(_score.index))
            .reset_index(drop=True))
    rnd['GeneID'] = _score['GeneID']
    _anat = anat.copy()
    rand_dir = os.path.join(anat['anat_dir'], str(seed))
    _anat.update({ 'anat_dir' : rand_dir,
        'outfile' : os.path.join(anat['anat_dir'], 'random')})
    _, go_rand = inference(exp, rnd, network_undirected, _anat, go)
    shutil.rmtree(rand_dir)
    return go_rand.dropna()

def inference_empiric_pvalue(exp, scores, network_undirected, anat, go, **kwargs):
    from joblib import Parallel, delayed
    perm = anat['permutations']
    subnetwork, go_scores = inference(exp, scores, network_undirected, anat, go)
    _anat = anat.copy()
    if 'viz_out' in _anat:
        del _anat['viz_out']
    MAX_INT = np.iinfo(np.int32).max
    seeds = np.random.randint(MAX_INT, size=perm)
    _scores = scores.reset_index(drop=True)
    rands = Parallel(n_jobs=35)(delayed(_rnd_go)(exp, _scores, network_undirected, _anat, go, seed) for seed in seeds)
    random = pd.concat(rands)[['term', 'p']].rename(columns={'p' : 'p_random'})
    df = (pd.merge(go_scores, random, how='left')
            .dropna(subset=['q'])
            .fillna(0))
    df['p>p_random'] = df['p'] > df['p_random']
    p_empiric = ((df.groupby('term')['p>p_random'].sum() / perm)
            .to_frame(name='p_empiric')
            .reset_index())
    result = pd.merge(go_scores, p_empiric)
    return subnetwork, result
