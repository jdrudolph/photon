import pandas as pd

def load(confidence, degree_threshold, db, ppi_network = None):
    """ read the network table from file
    filter for low confidence edges and high-degree nodes """
    ppi_network = db['ppi-network'] if ppi_network is None else ppi_network
    _anat = pd.read_table(ppi_network,
            names=['kin', 'sub', 'confidence', 'd'])
    anat = _anat[_anat['confidence'] > confidence]
    # Use degree cutoff to remove high-degree nodes
    _degk = anat.groupby('kin').apply(len)
    _degs = anat.groupby('sub').apply(len)
    _deg = pd.concat([_degk, _degs], axis=1).sum(axis=1)
    _low_deg = _deg[_deg < degree_threshold].index
    anat = anat[anat['kin'].isin(_low_deg) & anat['sub'].isin(_low_deg)].reset_index(drop=True)
    return anat

def to_directed(anat):
    """ duplicate undirected edges to form a directed network """
    found_in = ['found_in'] if 'found_in' in anat else []
    anat_directed = (pd.concat([anat, anat.rename(columns={'kin':'sub', 'sub':'kin'})])
            .reset_index(drop=True)
            [['kin', 'sub', 'confidence', 'd'] + found_in])
    anat_directed['d'] = 1
    return anat_directed

def add_unbiased_confidence(anat, evidence_path):
    """add large_scale evidence as confidence replacement to remove phosphositeplus bias
    
    >>> unbiased = add_unbiased_confidence(anat, 'db/anat/H_sapiens.txt')
    """
    evidence = pd.read_table(evidence_path)
    SOURCES = {'Y2H' : 'MI:0018'}
    evidence['found_in'] = evidence[list(SOURCES.values())].sum(axis=1)
    return pd.merge(anat, evidence, left_on=['kin', 'sub'],
            right_on=['interactor_a', 'interactor_b'])[['kin', 'sub', 'confidence', 'd', 'found_in']]


