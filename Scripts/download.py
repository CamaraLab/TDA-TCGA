import networkx
import multiprocessing
import json
import numpy
import numpy.linalg
import random
import requests
import scipy.stats
import pickle
import pylab


"""
GLOBAL METHODS
"""


def ParseAyasdiGraph(source, lab, user, password, name):
    """
    Parses Ayasdi graph given by the source ID and lab ID, and stores as name.gexf and name.pickle. user and
    password specify Ayasdi login credentials.
    """
    headers = {"Content-type": "application/json"}
    session = requests.Session()
    session.post('https://core.ayasdi.com/login', data={'username': user, 'passphrase': password})
    r = session.get('https://core.ayasdi.com/v0/sources/' + source + '/networks/' + lab)
    sp = json.loads(r.content)
    rows = [int(x['id']) for x in sp['nodes']]
    dic2 = {}
    for i in rows:
        payload = {"network_nodes_descriptions": [{"network_id": lab, "node_ids": [i]}]}
        r = session.post('https://core.ayasdi.com/v0/sources/' + source + '/retrieve_row_indices',
                         data=json.dumps(payload), headers=headers)
        dic2[i] = json.loads(r.content)['row_indices']
    with open(name + '.json', 'wb') as handle3:
        json.dump(dic2, handle3)
    rowcount = []
    with open(name + '.gexf', 'w') as g:
        g.write('<?xml version="1.0" encoding="UTF-8"?>\n')
        g.write('<gexf xmlns="http://www.gexf.net/1.2draft" version="1.2">\n')
        g.write('\t<graph mode="static" defaultedgetype="undirected">\n')
        g.write('\t\t<nodes>\n')
        for nod in sp['nodes']:
            g.write('\t\t\t<node id="' + str(nod['id']) + '" label="' + str(nod['row_count']) + '" />\n')
            rowcount.append(float(nod['row_count']))
        g.write('\t\t</nodes>\n')
        g.write('\t\t<edges>\n')
        for n5, edg in enumerate(sp['links']):
            g.write('\t\t\t<edge id="' + str(n5) + '" source="' + str(edg['from']) + '" target="' + str(edg['to'])
                    + '" />\n')
        g.write('\t\t</edges>\n')
        g.write('\t</graph>\n')
        g.write('</gexf>\n')

import argparse
parser = argparse.ArgumentParser()
parser.add_argument("source")
parser.add_argument("lab")
parser.add_argument("name")
args = parser.parse_args()

ParseAyasdiGraph(args.source,args.lab,'uer2102@columbia.edu','Columbia!',args.name)