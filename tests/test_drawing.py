import pytest
import phos
import phos.defaults
import phos.pipeline
import os
import pandas as pd
import phos.algo.subnetwork
from networkx.readwrite import json_graph
import json

@pytest.mark.quick
def test_drawing_small_example():
    exp = pd.DataFrame({
        "GeneID": [10, 20, 30, 40, 50],
        'Symbol': ['a','b','c','d','e'],
        'Amino.Acid': ['A', 'A', 'A', 'A', 'A'],
        'Position': [1, 2, 3, 4, 5],
        'avg': [1.0, 2, 3, 4, 5]})
    scores = pd.DataFrame({
        "GeneID": [1,2,3,4,5],
        "Significant": [True, True, True, False, True]})
    network = pd.DataFrame({"s": [10, 20, 30, 40], "t": [20, 30, 40, 50]})
    defaults = phos.defaults.make_defaults(os.path.abspath("."))
    G = phos.algo.subnetwork.create_graph(exp, scores, network, "test", defaults['db'], anchor=3)
    graph = json_graph.node_link_data(G)
    assert G.number_of_nodes() == len(graph['nodes']) == 5
    assert G.number_of_edges() == len(graph['links']) == 4
    assert type(graph['nodes']['id']) is str
    from jinja2 import Environment, FileSystemLoader
    template_dir = os.path.join(defaults['root'], 'templates')
    env = Environment(loader=FileSystemLoader(template_dir))
    graph_str = json.dumps(graph)
    HTML = env.get_template('result.html').render(task_id='test', graph=graph_str)
    with open('test.html', 'w') as f:
        f.write(HTML)
    assert False
