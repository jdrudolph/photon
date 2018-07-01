import pytest
import pandas as pd
from pandas.testing import assert_frame_equal
from phos.util.network import add_edge_is_directed

@pytest.mark.quick
def test_finding_edge_direction():
    network = pd.DataFrame({'kin': ['a', 'a', 'b', 'b', 'c', 'c'], 'sub': ['b', 'c', 'a', 'c', 'b', 'd']})
    expected = pd.DataFrame({'kin': ['a', 'a', 'b', 'c'], 'sub': ['b', 'c', 'c', 'd'], 'is_directed': [False, True, False, True]})
    actual = add_edge_is_directed(network).sort_values(by=['kin']).reset_index(drop=True)
    assert_frame_equal(expected, actual)
