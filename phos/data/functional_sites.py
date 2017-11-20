import pickle
import os

import pandas as pd

def filter_phospho(df):
    """ filter and format phosphosite extra tables """
    is_human = df['ORGANISM'] == 'human'
    ptms = df['MOD_RSD'].str.split('-')
    ptm_type = ptms.str.get(1)
    ptm_rsd = ptms.str.get(0)
    ptm_aa = ptm_rsd.str.get(0)
    ptm_pos = ptm_rsd.str.slice(1).astype(int)
    ptm_split = pd.concat([df['GENE_ID'], ptm_aa, ptm_pos], axis=1)
    phospho = ptm_split[(ptm_type == 'p') & is_human]
    phospho.columns = ['GeneID', 'Amino.Acid', 'Position']
    phospho = phospho.dropna()
    phospho['GeneID'] = phospho['GeneID'].astype(int)
    return phospho.drop_duplicates()

def _rename_amino_acid(df):
    """ 3 letter to single letter amino-acid code"""
    STY = {'Ser' : 'S', 'Thr' : 'T', 'Tyr' : 'Y'}
    amino_acid = STY.get(df['ORIGINAL'].iloc[0], float('NaN'))
    df['Amino.Acid'] = amino_acid
    return df

def _phosphosite_extra(filename, organism):
    """ parse phosphosite data """
    df = pd.read_csv(filename)
    if organism is not None:
        if hasattr(organism, '__iter__'):
            mask = df.ORGANISM.isin(organism)
        else:
            mask = df.ORGANISM == organism
        df.drop(df.index[~mask], inplace=True)
    return df

def regulatory_phosphosites(filename, organism=None, **kwds):
    return _phosphosite_extra(filename, organism)

def disease_associated_phosphosites(filename, organism=None, **kwds):
    return _phosphosite_extra(filename, organism)

def uniprot_humsavar(filename):
    """
    read amino-acid polymorphisms from
    http://uniprot.org/docs/humsavar
    """
    import re
    import io
    out = io.StringIO()
    out.write('\t'.join(['GENE_NAME', 'ACC', 'FTId', 'ORIGINAL',
        'POSITION', 'MUTATED', 'VARIANT', 'dbSNP', 'DISEASE']) + '\n')
    with open(filename) as f:
        s = re.compile('(\S+)\s+(\S+)\s+(\S+)\s+p\.([a-zA-Z]+)(\d+)([a-zA-Z]+)\s+(\S+)\s+(\S+)\s+(.*)') # one regex to rule them all
        for i,line in enumerate(f):
            if i < 30: # header
                continue
            if line == '\n': # reached end
                break
            out.write('\t'.join(s.match(line.strip()).groups()) + '\n')
    out.seek(0)
    return pd.read_table(out, na_values='-')

def read_files(db):
    """ return functional sites sources
    dis - disease associated sites from uniprot
    reg - regulatory sites from uniprot
    muts - SAPs
    """
    import phos.data.mapping as mapping
    dis = filter_phospho(disease_associated_phosphosites(db['phosphosite']['disease_associated_sites']))
    reg = filter_phospho(regulatory_phosphosites(db['phosphosite']['regulatory_sites']))
    _muts = (uniprot_humsavar(db['uniprot_humsavar']).dropna()
            .groupby('ORIGINAL').apply(_rename_amino_acid).dropna()
            .rename(columns={'POSITION': 'Position'}))
    muts = mapping.map_protein_groups(_muts, 'ACC')[['GeneID', 'Amino.Acid', 'Position']]
    return dis, reg, muts

def functional_sites(db):
    _cache = '.functional_sites.pkl'
    if os.path.isfile(_cache):
        with open(_cache, 'rb') as f:
            functional_sites = pickle.load(f)
    else:
        dis, reg, muts = read_files()
        functional_sites = (dis.append(reg).append(muts).drop_duplicates())
        with open(_cache, 'wb') as f:
            pickle.dump(functional_sites, f)
