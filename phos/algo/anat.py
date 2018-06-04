import os
import subprocess

import pandas as pd
import networkx as nx
from networkx.readwrite import json_graph
import json

def run_locally(network_undirected, terminals, data, anchor,
        alpha, task_id, work_dir, viz_out=None, steinprt='steinprt', **kwargs):
    """
    run anat subnetwork inference ** requires steinprt binary **
    
    :param network_undirected: network
    :param terminals:
    :param dataset: pandas.DataFrame with columns: Symbol, GeneID, Amino.Acid, Position, avg
    :param anchor: anchor
    :param alpha: anat parameter, tradeoff between shortest path and steiner tree
    :param viz_out: file path. Result will be saved at `viz_out`.js and `viz_out`.html
    :returns subnetwork: the resulting subnetwork as a table
    """
    anat_dir = os.path.join(work_dir, task_id, 'anat')
    os.makedirs(anat_dir, exist_ok=True)
    _terminals = terminals.to_frame(name='terminal').reset_index(drop=True)
    _terminals['anchor'] = anchor
    _terminals = _terminals[_terminals['terminal'] != anchor][['anchor', 'terminal']]
    (network_undirected.iloc[:,:4]
        .to_csv(os.path.join(anat_dir, 'network'), index=False, sep='\t', header=None))
    _terminals.to_csv(os.path.join(anat_dir, 'terminals'), index=False, sep='\t', header=None)
    outfile = os.path.join(anat_dir, 'anat.out')
    output = subprocess.check_output([steinprt,
        '-n', 'network',
        '-s', 'terminals',
        '-b', '{:0.2f}'.format(alpha),
        '-f', anat_dir,
        '-r', outfile]) # -f flag ignored by steinprt
    with open(outfile + '.log', 'wb') as f:
        f.write(output)
    subnetwork = pd.read_table(outfile, sep=' ',
            names=['s', 't', 'confidence', 'transformed', 'directed'])
    return subnetwork


"""
                  CALLING ANAT SERVER FROM PYTHON EXAMPLE

    import uuid
    sessionId = str(uuid.uuid4())
    anchor = 1950
    terminals = [1120, 112, 1956]

    from anat import network
    nodes, edges = network(sessionId, anchor, terminals)

"""
import requests
import time
import xml.etree.ElementTree as ET

def get_result(sessionID):
    headers = {'Accept' : 'text/xml, multipart/related',
            'Content-Type' : 'text/xml; charset=utf-8',
            'SOAPAction' : '"getResult"'}
    data = """<?xml version='1.0' encoding='UTF-8'?>
    <S:Envelope xmlns:S="http://schemas.xmlsoap.org/soap/envelope/">
        <S:Body><ns2:sessionId xmlns:ns2="network">{}</ns2:sessionId></S:Body>
    </S:Envelope>""".format(sessionID)
    result_response = requests.post('http://anat.cs.tau.ac.il/AnatWeb/AnatServer', headers=headers, data=data)
    if not result_response.ok:
        raise ValueError('Response status not OK',
                '\nRequest', result_response.request.body,
                '\nResponse', result_response.content)
    return result_response.text

def add_edges(network_undirected):
    edgesData = ["""       <ns2:edgesData>
         <ns2:action>SET_UNDIRECTED</ns2:action>
         <ns2:additionalInfo>edge</ns2:additionalInfo>
         <ns2:confidence>{confidence}</ns2:confidence>
         <ns2:fromNodeId>{source}</ns2:fromNodeId>
         <ns2:toNodeId>{target}</ns2:toNodeId>
       </ns2:edgesData>""".format(source=source,
           target=target, confidence=confidence)
       for (source, target), confidence in
       zip(network_undirected[['kin', 'sub']].values,
           network_undirected['confidence'])]
    return '\n'.join(edgesData)
    
def make_terminals(nodes):
    return '\n'.join([
        '      <ns2:terminals>{}</ns2:terminals>'.format(node)
        for node in nodes])

def make_set(nodes):
    return '\n'.join([
        '      <ns2:set>{}</ns2:set>'.format(node) for node in nodes])

def submit_job_no_anchor(sessionId, network_undirected, terminals):
    edgesData = add_edges(network_undirected)
    setData = make_set(terminals)
    headers = { 'SOAPAction' : '"calculateProjectionSubNetwork"',
            'charset' : 'utf-8',
            'content-type' : 'text/xml' }
    data = """<?xml version="1.0" encoding="UTF-8"?>
    <S:Envelope xmlns:S="http://schemas.xmlsoap.org/soap/envelope/">
      <S:Body>
        <ns2:projectionParams xmlns:ns2="network">
          <ns2:algorithmType>PROJECTIONANALYSIS</ns2:algorithmType>
          <ns2:backGroundNetwork>
            <ns2:networkName>my_network</ns2:networkName>
            {edgesData}
          </ns2:backGroundNetwork>
          <ns2:baseNetworkFileName>E_empty.net</ns2:baseNetworkFileName>
          <ns2:edgeConstraintMap/>
          <ns2:nodeConstraintMap/>
          <ns2:sessionId>{sessionId}</ns2:sessionId>
          <ns2:title>mytitle</ns2:title>
          <ns2:lengthPenalty>25</ns2:lengthPenalty>
          <ns2:margin>0</ns2:margin>
          <ns2:algorithmType>PROJECTIONANALYSIS</ns2:algorithmType>
          <ns2:granularity>100</ns2:granularity>
          {setData}
          <ns2:subAlgorithm>CLUSTERING</ns2:subAlgorithm>
        </ns2:projectionParams>
      </S:Body>
    </S:Envelope>
    """.format(edgesData=edgesData, setData=setData, sessionId=sessionId)
    submit_response = requests.post(
            'http://anat.cs.tau.ac.il/AnatWeb/AnatServer',
            headers=headers, data=data)
    if not submit_response.ok:
        raise ValueError('Response status not OK',
                '\nRequest', submit_response.request.body,
                '\nResponse', submit_response.content)
    return submit_response.text 

def submit_job(sessionId, network_undirected, anchor, terminals):
    edgesData = add_edges(network_undirected)
    terminalsData = make_terminals(terminals)
    headers = { 'SOAPAction' : '"calculateExplanatorySubNetwork"',
            'charset' : 'utf-8',
            'content-type' : 'text/xml' }
    data = """<?xml version='1.0' encoding='UTF-8'?>
      <S:Envelope xmlns:S="http://schemas.xmlsoap.org/soap/envelope/">
        <S:Body>
        <ns2:explanatoryParameters xmlns:ns2="network">
          <ns2:algorithmType>EXPLANATORYPATHWAYS</ns2:algorithmType>
          <ns2:backGroundNetwork>
            <ns2:networkName>my_network</ns2:networkName>
            {edgesData}
          </ns2:backGroundNetwork>
          <ns2:baseNetworkFileName>E_empty.net</ns2:baseNetworkFileName>
          <ns2:edgeConstraintMap/>
          <ns2:nodeConstraintMap/>
          <ns2:sessionId>{sessionId}</ns2:sessionId>
          <ns2:title>my_title</ns2:title>
          <ns2:lengthPenalty>25</ns2:lengthPenalty>
          <ns2:margin>0</ns2:margin>
          <ns2:algorithmType>EXPLANATORYPATHWAYS</ns2:algorithmType>
          <ns2:alpha>0.25</ns2:alpha>
          <ns2:anchors>{anchor}</ns2:anchors>
          <ns2:predictTF>false</ns2:predictTF>
          <ns2:subAlgorithm>APPROXIMATE</ns2:subAlgorithm>
          {terminalsData}
          <ns2:terminalsToAnchors>false</ns2:terminalsToAnchors>
        </ns2:explanatoryParameters>
        </S:Body>
      </S:Envelope>
    """.format(edgesData=edgesData, sessionId=sessionId,
            anchor=anchor, terminalsData=terminalsData)
    submit_response = requests.post(
            'http://anat.cs.tau.ac.il/AnatWeb/AnatServer',
            headers=headers, data=data)
    if not submit_response.ok:
        raise ValueError('Response status not OK',
                '\nRequest', submit_response.request.body,
                '\nResponse', submit_response.content)
    return submit_response.text 

def parse_result(result):
    """ Parse ANAT server response. Returns pd.DataFrame with edges or None if no edges were found
    """
    network_graph = result[0][0]
    edges = {}
    for x in network_graph:
        if x.tag == '{network}edges':
            [edges.setdefault(a.tag.replace('{network}', ''), [])
                    .append(a.text) for a in x]
    edges = pd.DataFrame(edges)
    try:
        subnetwork = (edges[['id1', 'id2']]
            .rename(columns={'id1' : 's', 'id2' : 't'}))
        return subnetwork
    except KeyError:
        return None


def network(sessionId, terminals, anchor=None, **kwargs):
    if anchor is None:
        submit_response = submit_job_no_anchor(sessionId, terminals, **kwargs)
    else:
        submit_response = submit_job(sessionId, anchor, terminals, **kwargs)
    got_result = False
    while not got_result:
        result = ET.fromstring(get_result(sessionId))
        got_result = len(result[0][0]) > 0
        time.sleep(0.5)
    return parse_result(result)


def remote_network(sessionId, network_undirected, terminals, anchor=None, **kwargs):
    if anchor is None:
        submit_response = submit_job_no_anchor(sessionId, network_undirected, terminals, **kwargs)
    else:
        submit_response = submit_job(sessionId, network_undirected, anchor, terminals, **kwargs)
    retries = 500 # more than half an hour
    got_result = False
    while not got_result and retries > 0:
        result = ET.fromstring(get_result(sessionId))
        got_result = len(result[0][0]) > 0
        time.sleep(5)
        retries = retries - 1
    if retries <= 0:
        print("Failed to obtain ANAT results. Exceeded retry limit of 500", flush=True)
        return None
    return parse_result(result)

