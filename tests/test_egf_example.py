import uuid
import pytest
import phos
import phos.defaults
import phos.pipeline
import os
import pandas as pd

parameters = {
        "activity" : {
            "min_size" : 4,
            "permutations" : 1000,
            "side" : "greater"
            },
        "anat" : {
            "alpha" : 0.25,
            "anchor" : -1
            },
        "go" : {
            "max_category_size" : 500,
            "min_category_size" : 5,
            "max_category_depth" : 10
            },
        "ppi-network" : {
            "confidence" : 0.5,
            "degree_threshold" : 150
            }
        }

@pytest.mark.parametrize('clean_dir_with_data', ['static/data.csv'], indirect=['clean_dir_with_data'])
def test_egf_data_example(clean_dir_with_data):
    defaults = phos.defaults.make_defaults(os.path.abspath('.'))
    _parameters = parameters.copy()
    _parameters['anat']['anchor'] = 1950
    results = phos.pipeline._run(str(uuid.uuid4()), 'data.csv', _parameters, defaults['db']) 

@pytest.mark.parametrize('clean_dir_with_data', ['static/data.csv'], indirect=['clean_dir_with_data'])
def test_one_site_per_gene(clean_dir_with_data):
    defaults = phos.defaults.make_defaults(os.path.abspath('.'))
    _parameters = parameters.copy()
    _parameters['anat']['anchor'] = 1950
    df = pd.read_csv('data.csv')
    df.drop_duplicates(subset=['GeneID']).to_csv('one_site_per_gene.csv', index=False)
    results = phos.pipeline._run(str(uuid.uuid4()), 'one_site_per_gene.csv', _parameters, defaults['db'])

@pytest.mark.parametrize('clean_dir_with_data', ['static/data.csv'], indirect=['clean_dir_with_data'])
def test_tiny_dataset(clean_dir_with_data):
    defaults = phos.defaults.make_defaults(os.path.abspath('.'))
    _parameters = parameters.copy()
    _parameters['anat']['anchor'] = 1950
    df = pd.read_csv('data.csv')
    df.head(100).to_csv('tiny_dataset.csv', index=False)
    with pytest.raises(ValueError) as exc:
        results = phos.pipeline._run(str(uuid.uuid4()), 'tiny_dataset.csv', _parameters, defaults['db'])
        assert 'No significant signaling functionality scores found' in str(exc.value)
