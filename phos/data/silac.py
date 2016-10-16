""" utility function for SILAC (log-fold) based datasets """
import pandas as pd
import numpy as np
from statsmodels.stats.multitest import fdrcorrection
from scipy.stats import ttest_1samp
import phos.util.mapping as mapping

def threshold(data, fold_changes, t=2, **kwargs):
    return pd.DataFrame(
            (data[fold_changes] > t).all(1) | (data[fold_changes] < -t).all(1),
            columns=['rej'])

def one_sample_ttest(grouped, data, fold_changes, **kwargs):
    """ one sample ttest for differential phosphorylation of a gene,
    fdr correction at alpha 0.05 """
    data_mean = np.mean(data[fold_changes].values)
    test_results = (grouped
            .apply(lambda df: ttest_1samp(np.hstack(df[fold_changes].values), data_mean)))
    test_T = test_results.str.get(0)
    test_p = test_results.str.get(1)
    rej, q = fdrcorrection(test_p.values)
    df = pd.concat([test_T, test_p], axis=1)
    df.columns = ['T', 'p']
    df['rej'] = rej
    df['q'] = q
    return df.reset_index()

def experiment(df, fold_changes, uniprot_mapping, geneinfo, **kwargs):
    """ use 'experiments' entry to map df to GeneID and extract relevant columns
   
    >>> mapped = silac.experiment(df, **params)

    :param df: pandas.DataFrame to extract from
    :param fold_changes: list of data columns
    :param uniprot_mapping: string specifying the uniprot mapping file location
    :param geneinfo: string specifying the geninfo file location
    :returns mapped: pandas.DataFrame with the experiment mapped to geneids
    """
    site_cols = ['Amino.Acid', 'Position']
    df = (df.dropna(subset=site_cols)
            .dropna(subset=fold_changes, how='all'))
    df = df.rename(columns={fc : 'fc' for fc in fold_changes})
    geneid2symbol = pd.read_table(geneinfo, compression='gzip', comment='#',
            usecols = [1, 2], names=['GeneID', 'Symbol'])
    mapped = (mapping
            .map_protein_groups(df, 'Uniprot', 'GeneID', uniprot_mapping)
            [['GeneID'] + site_cols + ['fc']]
            .dropna()
            .drop_duplicates(subset = ['GeneID'] + site_cols)
            .merge(geneid2symbol))
    mapped['Position'] = mapped['Position'].astype(int)
    mapped['fc'] = (mapped['fc'] - mapped['fc'].mean()) / mapped['fc'].std()
    if type(mapped['fc']) is pd.Series:
        mapped['avg'] = mapped['fc']
    else: # type is pd.DataFrame
        mapped['avg'] = mapped['fc'].mean(axis=1)
    return mapped
