""" utility functions for id mapping """
import os.path

import pandas as pd

def map_protein_groups(data, uniprot, geneid, uniprot_mapping, astype=int):
    """ turn protein group dataset into gene datset

    >>> genes = map_protein_groups(data, 'Uniprot ID', 'GeneID', 'idmapping.dat.gz')
    
    :param data: pandas.DataFrame to be mapped
    :param uniprot: string with the uniprot column name
    :param geneid: string specifying the target geneid column name,
                    has to correspond to idmapping.dat db name
    :param uniprot_mapping: string with the path to the uniprot id-mapping file
    :param astype: convert mapped ids to type
    :returns data: mapped"""
    uni2gid = dict(load_uniprot(uniprot_mapping, geneid).values)
    data[geneid] = data[uniprot].apply(_mapper, args=(uni2gid,))
    data = data.dropna(subset=[geneid])
    data = split_and_stack(data, geneid, ';')
    data[geneid] = data[geneid].astype(astype)
    return data

def _mapper(unis, uni2gid):
    """ maps uniprot to geneid and joins on ';' """
    if type(unis) is float:
        return uni2gid.get(unis, None)
    mapped = {uni2gid.get(u, None) for u in unis.split(';')}
    filtered = [m for m in mapped if m is not None]
    if len(filtered) == 0:
        return None
    return ';'.join(filtered)



def load_uniprot(filename, db='GeneID', force=False):
    """
    read unprot mapping from file. cache uniprot table to save time
    
    :param filename: gzipped uniprot mapping file path HUMAN_9606_idmapping.dat.gz
    :param force: force reparsing the file
    :returns uniprot: pandas.DataFrame mapping table
    """
    cache = '{}_{}_.pkl'.format(filename, db)
    if os.path.isfile(cache) and not force:
        uniprot = pd.read_pickle(cache)
    else: 
        uniprot = parse_uniprot_mapping(filename, db)
        uniprot.to_pickle(cache)
    return uniprot

def parse_uniprot_mapping(filename, db='GeneID'):
    """
    read uniprot mapping from file
    
    :param filename: gzipped uniprot mapping file path HUMAN_9606_idmapping.dat.gz
    """
    raw_data = pd.read_table(filename, compression='gzip',
            names=['Uniprot', 'db', 'dbid'])
    geneid_data = raw_data[raw_data['db'] == db]
    mapping = pd.pivot_table(geneid_data, 'dbid', index='Uniprot', columns='db',
            aggfunc=lambda x : ';'.join(x)).reset_index()
    return mapping

def split_and_stack(frame, colname, sep=None):
    """
    Converts a row of delimited values in a dataframe to multiple rows.
    The resulting table is useful for joins.

    >>> import pandas as pd 
    >>> df = pd.DataFrame.from_items([('col1', ['1 2 3']), ('col2', ['a'])])
    >>> splitted = split_and_stack(df, 'col1')
    >>> splitted.shape
    (3, 2)

    >>> splitted
      col2 col1
    0    a    1
    0    a    2
    0    a    3

    :param frame: the data frame to process
    :param colname: the name of the column containing the delimited data
    :param sep: an optional pattern for split(). if not specified split on whitespace
    """
    mask = frame[colname].isnull()

    masked_frame=frame[~mask]
    mapping = pd.DataFrame([{i: v for i,v in enumerate(x)} for x in
                             masked_frame[colname].str.split(sep).tolist()],
                             index=masked_frame.index).stack().reset_index(1,drop=True)
    mapping.name = colname # required for the merge

    return frame.drop(colname, axis=1).join(mapping)

