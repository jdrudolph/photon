import pandas as pd
import numpy as np

from statsmodels.stats.multitest import fdrcorrection

def preprocess(exp, network, min_size):
    _kins = (pd.merge(exp, network, left_on='GeneID', right_on='sub')
            .groupby('kin')
            .apply(len))
    kins = _kins[_kins > min_size]
    filtered_network = network[network['kin'].isin(kins.index)]
    return kins, filtered_network

def run_one_sample_regular_weights(exp, network, beta=0):
    """
    weighted mean and standard deviation
    inspired by http://www.gnu.org/software/gsl/

    FIXME: remove beta
    """
    df = (pd.merge(exp, network, left_on='GeneID', right_on='sub')
            .rename(columns={'avg' : 'x', 'confidence' : 'w'}))
    if df.shape[0] == 0:
	    raise ValueError("None of the 'GeneID' could be matched to the network. Only human entrez GeneIDs are supported.")
    df['wx'] = df['w'] * df['x']
    df['w(x^2)'] = df['w'] * np.power(df['x'], 2)
    df['w^2'] = np.power(df['w'], 2)
    kin = df.groupby('kin')
    n = kin.apply(len)
    V1 = kin['w'].sum()
    V2 = kin['w^2'].sum()
    var = kin['w(x^2)'].sum() / (V1 - (V2 / V1))
    std = np.sqrt(var)
    mean = kin['wx'].sum() / V1
    t = mean / (np.power(std, beta) / np.sqrt(n))
    return t.to_frame(name='t')

def empiric(exp, network, min_size, permutations, side):
    num_neighbors, fnetwork = preprocess(exp, network, min_size)
    run = run_one_sample_regular_weights
    _exp = exp.reset_index(drop=True)
    emp = run(_exp, fnetwork)
    emp['pos'] = 0
    emp['neg'] = 0
    for _ in range(permutations):
        rnd = (_exp['avg']
                .reindex(np.random.permutation(_exp.index))
                .reset_index(drop=True))
        if type(rnd) is pd.Series:
            rnd = rnd.to_frame(name='avg')
        rnd['GeneID'] = _exp['GeneID']
        rnd_act = run(rnd, fnetwork)
        emp['pos'] = emp['pos'] + (emp['t'] < rnd_act['t'])
        emp['neg'] = emp['neg'] + (emp['t'] >= rnd_act['t'])
    emp['p_greater'] = emp['pos'] / permutations
    emp['p_lesser'] = emp['neg'] / permutations
    emp['p_twosided'] = emp[['p_greater', 'p_lesser']].min(axis=1)
    rej, q = fdrcorrection(emp['p_greater'].append(emp['p_lesser']))
    n = emp.shape[0]
    emp['rej_greater'] = rej[:n]
    emp['rej_lesser'] = rej[n:]
    emp['rej_twosided'] = rej[:n] | rej[n:]
    
    emp['q_greater'] = q[:n]
    emp['q_lesser'] = q[n:]
    emp['q_twosided'] = np.min([q[:n], q[n:]], axis=0)

    mask = emp['p_greater'] < emp['p_lesser']
    # continuity correction with 1/permutations to avoid log(0) = -Inf
    emp['score_empiric'] = ((mask * -np.log(emp['p_twosided'] + 1/permutations))
            + (~mask * np.log(emp['p_twosided'] + 1/permutations)))
    if side not in {'greater', 'twosided', 'lesser'}:
        raise ValueError('side parameters has to be one of "greater", "twosided", "less", was {}'.format(side))
    emp['Significant'] = emp['rej_{}'.format(side)]
    return (emp.drop(['pos', 'neg'], 1)
            .reset_index()
            .rename(columns={'kin' : 'GeneID'}))
