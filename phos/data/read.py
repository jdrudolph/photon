"""
tools for reading configuration info and data
"""
import os

import pandas as pd

import phos.data
import phos.util.io as io

def load(input, setting, experiment, uniprot_mapping, geneinfo, **kwargs):
    """
    read the dataset from file as specified

    >>> dataset, aux = config('egf.py')
    >>> exp = load(dataset, aux)

    :param dataset: dictionary
    :param aux: dictionary
    :returns exp: pandas.DataFrame
    """
    cache = experiment.get('cache', None)
    if cache is not None and os.path.isfile(cache):
        exp = pd.read_pickle(cache)
    else:
        df = io.read(**input)
        experiment_func = phos.data.get_experiment(setting)
        exp = experiment_func(df, uniprot_mapping=uniprot_mapping, geneinfo=geneinfo,
                **experiment)
        if cache is not None:
            exp.to_pickle(cache)
    return exp

def config(configfile):
    """
    read a configfile and return the dataframe
    
    >>> configs = config('egf.py')

    :param configfile: the path to the configuration file
    :returns dataset: the dictionary with the configu info
    """
    configs = {}
    with open(configfile) as f:
        exec(f.read(), configs)
    return configs
