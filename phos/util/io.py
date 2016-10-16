""" reading and writing data """
import pandas as pd
import os

def read(path, format, NaN=None, data=None, query=None, cache=None, **kwargs):
    """ read dataset raw file
    >>> df = read(**dataset['file']) """
    if cache is not None:
        if os.path.isfile(cache):
            raw = pd.read_pickle(cache)
        else:
            raw = _read(path, format, NaN, data, query, **kwargs)
            raw.to_pickle(cache)
    else:
        raw = _read(path, format, NaN, data, query, **kwargs)
    return raw

def _read(path, format, NaN=None, data=None, query=None, **kwargs):
    if format == 'tsv':
        raw = pd.read_table(path, **kwargs)
    elif format == 'xlsx':
        raw = pd.read_excel(path, **kwargs)
    else:
        raise NotImplementedError('unknown data format', format)
    if NaN is not None:
        raw = raw.mask(raw == NaN)
    if data is not None:
        raw = raw.dropna(subset=data)
    if query is not None:
        raw = raw.query(query)
    return raw
