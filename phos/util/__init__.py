import phos.util.mapping
import phos.util.io
import json
from collections import OrderedDict

def read_json(filename):
    with open(filename) as f:
        return json.load(f, object_pairs_hook=OrderedDict)

def deep_update(old, new):
    """ updates the second level of a dictionary

    >>> d = {'a' : {'b' : 1, 'c': 2}}
    >>> deep_update(d, {'a' : {'b' : 2}})
    >>> d
    {'a' : {'b' : 2, 'c' : 2}}
    """
    for k, v in new.items():
        old[k].update(v)

def num(string):
    """ convert string to int/float
    
    >>> type(num(1))
    int
    >>> type(num(1.5))
    float
    """
    try:
        return int(string)
    except ValueError:
        return float(string)

def uuid():
    import uuid
    return str(uuid.uuid4())
