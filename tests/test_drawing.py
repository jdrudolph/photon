import pytest
import phos
import phos.defaults
import phos.pipeline
import os
import pandas as pd
import phos.algo.subnetwork

@pytest.mark.quick
def test_drawing_small_example():
    exp = pd.DataFrame({"GeneID": [1,2,3,4,5]})
    scores = pd.DataFrame({"GeneID": [1,2,3,4,5], "Significant": [True, True, True, False, True]})
    network = pd.DataFrame({"s": [1,2,3,4], "t": [2,3,4,5]})
    defaults = phos.defaults.make_defaults(os.path.abspath("../."))
    G = phos.algo.subnetwork.create_graph(exp, scores, network, "test", defaults['db'], anchor=3)
