import uuid
import pytest
import phos
import phos.defaults
import phos.pipeline
import os
parameters = {
        "activity" : {
            "min_size" : 4,
            "permutations" : 100
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
