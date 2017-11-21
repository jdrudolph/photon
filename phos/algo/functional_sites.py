import os

import pandas as pd

import phos.algo.logistic_regression as clf
import phos.util.mapping as mapping

def predict_functional_sites(exp, scores, db, _functional_features = ['num_sites'], **kwargs):
    num_sites = get_sites_per_geneid(db['uniprot_mapping'],
            **db['phosphosite'])
    predicted_functional = (exp.merge(scores, on='GeneID')
            .merge(num_sites, on='GeneID'))
    known_functional = get_functional_sites(**db['phosphosite'])
    known_functional['functional'] = 1.0
    merged = pd.merge(predicted_functional, known_functional, how='left')
    predicted_functional['label'] = merged['functional'].notnull()
    features = ['avg'] + _functional_features + ['score_empiric']
    if merged.shape[0] > 10:
        res = clf.logistic_regression(predicted_functional, features)
        probas, coefs, mean_fpr, mean_tpr, mean_auc, classifier = res
    else:
        probas, coefs, mean_fpr, mean_tpr, mean_auc, classifier = (float('NaN'), [], float('NaN'), float('NaN'), float('NaN'), None)
    predicted_functional['proba'] = probas
    return predicted_functional, coefs, mean_fpr, mean_tpr, mean_auc, classifier

def get_sites_per_geneid(uniprot_mapping, all_sites, **kwargs):
    """ calculate how many known sites are on each geneid.
    based on all sites reported in phosphositeplus

    :param all_sites: path to the phosphorylation site dataset
    :param uniprot_mapping: path to uniprot mapping file
    :returns num_sites:
    """
    cache = all_sites + '.pkl'
    if os.path.isfile(cache):
        mapped_sites = pd.read_pickle(cache)
    else:
        sites = pd.read_table(all_sites, compression='gzip', skiprows=3)
        human_sites = sites[sites['ORGANISM'] == 'human']
        mapped_sites = mapping.map_protein_groups(human_sites, 'ACC_ID',
                'GeneID', uniprot_mapping)
        mapped_sites.to_pickle(cache)
    num_sites = (mapped_sites[['GeneID', 'MOD_RSD']]
            .drop_duplicates()
            .groupby('GeneID')
            .apply(len)
            .to_frame(name='num_sites')
            .reset_index())
    return num_sites

def get_functional_sites(regulatory_sites, disease_associated_sites, **kwargs):
    """ load known functional sites from phosphosite
    (i) regulatory sites (ii) disease-associated sites

    :param regulatory_sites: file path
    :param disease_associated_sites: file path
    """
    reg = _format_phosphosite_table(pd.read_csv(regulatory_sites))
    da = _format_phosphosite_table(pd.read_csv(disease_associated_sites))
    functional = pd.concat([reg, da]).drop_duplicates()
    return functional

# HELPER FUNCTIONS

def _format_phosphosite_table(df):
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
