import pytest
import phos
import phos.defaults
import phos.pipeline

@pytest.mark.parametrize('clean_dir_with_data', ['static/data.csv'], indirect=['clean_dir_with_data'])
def test_egf_data_example(clean_dir_with_data):
    parameters = phos.defaults.parameters.copy()
    phos.defaults.db = phos.defaults.make_db('.')
    parameters['anat']['anchor'] = 1950
    phos.pipeline.run("test", "data.csv", parameters) 
