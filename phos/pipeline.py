#!/usr/bin/env python
import os
from os.path import join
import pandas as pd
import numpy as np
import json

import phos.util

import phos.data.read as reader
import phos.data.network as networker

import phos.algo.activity as activity
import phos.algo.subnetwork as subnetwork
import phos.algo.functional_sites as functional_sites

from phos.defaults import make_defaults

def print_progress(message):
    print('progress:', message)

def _run(task_id, data, parameters, db, set_progress=print_progress):
    exp = pd.read_csv(data)
    columns = ['GeneID', 'Amino.Acid', 'Position', 'avg', 'Symbol']
    if len(exp.columns) < len(columns):
	    raise ValueError('Found {} columns but expected {}. Please check that the input format is correct!'.format(len(exp.columns), len(columns)))
    if not (exp.columns == columns).all():
        raise ValueError('Column names are not {}'.format(columns))
    if not exp['GeneID'].dtypes == np.int:
        raise ValueError('GeneID could not be parsed as integer, make sure the column does not contain "NaN" or ";"')
    set_progress('1/5 read data file, constructing network...')
    network_undirected = networker.load(db=db, **parameters['ppi-network'])
    network = networker.to_directed(network_undirected)
    set_progress('2/5 constructed the network, calculating functionality scores...')
    scores = activity.empiric(exp, network, **parameters['activity'])
    if len(scores[scores['Significant']]['GeneID']) < 1:
        raise ValueError('No significant signaling functionality scores found, possibly due to insufficient data.')
    set_progress('3/5 calculated scores, infering subnetwork using ANAT...')
    subnet, go_scores = subnetwork.inference(exp, scores,
            network_undirected, db=db, task_id=task_id, **parameters)
    set_progress('4/5 inferred subnetwork, making functional site predictions')
    predictions = functional_sites.predict_functional_sites(exp, scores, db=db,
            **parameters)
    set_progress('5/5 all done')
    return exp, scores, subnet, go_scores, predictions

def run(task_id, data, parameters, template_dir, work_dir, db, set_progress=print_progress):
    if parameters['anat'].get('anchor', None) < 1:
        del parameters['anat']['anchor']
    exp, scores, subnet, go_scores, predictions = _run(task_id, data, parameters, db=db, set_progress=set_progress)
    os.makedirs(join(work_dir, task_id), exist_ok=True)
    with open(join(work_dir, task_id, 'parameters.json'), 'w') as f:
        json.dump(parameters, f, indent=1)
    scores.to_csv(join(work_dir, task_id, 'scores.csv'), index=False)
    subnet.to_csv(join(work_dir, task_id, 'subnet.csv'), index=False)
    go_scores.to_csv(join(work_dir, task_id, 'go_scores.csv'), index=False)
    HTML = subnetwork.draw(exp, scores, subnet, task_id, template_dir=template_dir,
            anchor=parameters['anat'].get('anchor', None), db=db)
    with open(join(work_dir, task_id, 'result.html'), 'w') as f:
        f.write(HTML)
    predictions[0].to_csv(join(work_dir, task_id, 'predictions.csv'), index=False)

if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(description='PHOTON: PHOsphoproteomic dissecTiOn using Networks')
    parser.add_argument('data',
            help='data file')
    parser.add_argument('task_id', type=str, default=phos.util.uuid(),
            help='determines the save location, will be generated')
    parser.add_argument('parameters', help='read parameters from json')
    parser.add_argument('-r', '--root-dir', default=os.getcwd(), help='root directory')
    parser.add_argument('-t', '--template-dir', default=os.path.join(os.getcwd(), 'templates'), help='templates directory')
    parser.add_argument('-a', '--anchor', type=int, default=-1,
            help='GeneID for the subnetwork reconstruction anchor')
    parser.add_argument('-p', '--parameter', action='append',
            default = [],
            help='change parameter --parameter alpha=0.1, \
                    can be supplied mutliple times')
    args = parser.parse_args()
    defaults = make_defaults(args.root_dir)
    parameters = phos.util.read_json(args.parameters)
    parameters['anat']['anchor'] = args.anchor if args.anchor > 0 else None
    # update default values
    for _param in args.parameter:
        k,_v = _param.split('=')
        v = phos.util.num(_v)
        for cat in parameters:
            if k in parameters[cat]:
                parameters[cat][k] = v
                print(cat, k, v)
                break
    task_id = args.task_id
    os.makedirs(os.path.join(defaults['work'], task_id), exist_ok=True)
    print('running', task_id)
    run(task_id, args.data, parameters, work_dir=defaults['work'], template_dir=args.template_dir, db=defaults['db'])
