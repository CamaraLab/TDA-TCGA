import networkx
import multiprocessing
import json
import numpy
import numpy.linalg
import random
import requests
import pickle

def ParseAyasdiGraph(source, network, user, password, name):
    """
    Parses Ayasdi graph given by the source ID and network ID, and stores as name.gexf and name.pickle. user and
    password specify Ayasdi login credentials.
    """
    r = session.get('https://core.ayasdi.com/v1/sources/' + source + '/networks/' + network)
    sp = json.loads(r.content)
    rows = [int(x['id']) for x in sp['nodes']]
    dic2 = {}
    for i in rows:
        payload = {"network_nodes_descriptions": [{"network_id": network, "node_ids": [i]}]}
        r = session.post('https://core.ayasdi.com/v1/sources/' + source + '/retrieve_row_indices',
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
            g.write('\t\t\t<node id="' + str(nod['id']) + '" networkel="' + str(nod['row_count']) + '" />\n')
            rowcount.append(float(nod['row_count']))
        g.write('\t\t</nodes>\n')
        g.write('\t\t<edges>\n')
        for n5, edg in enumerate(sp['links']):
            g.write('\t\t\t<edge id="' + str(n5) + '" source="' + str(edg['from']) + '" target="' + str(edg['to'])
                    + '" />\n')
        g.write('\t\t</edges>\n')
        g.write('\t</graph>\n')
        g.write('</gexf>\n')


headers = {"Content-type": "application/json"}
session = requests.Session()
session.post('https://core.ayasdi.com/login', data={'username': 'uer2102@columbia.edu', 'passphrase': 'ColumbiaAyasdi2015!'})


import argparse
parser = argparse.ArgumentParser()
parser.add_argument("cancer_name")
parser.add_argument("source")
parser.add_argument("column_set_id")
parser.add_argument("metric")
parser.add_argument("lens")
parser.add_argument("res_start")
parser.add_argument("res_stop")
parser.add_argument("res_interval")
parser.add_argument("gain_start")
parser.add_argument("gain_stop")
parser.add_argument("gain_interval")
args = parser.parse_args()

if args.metric=="Cor" :
	metric="Correlation"

if args.lens=="Neigh" :
	len1="Neighborhood Lens 1"
	len2="Neighborhood Lens 2"



dict={}
gain_range=list(numpy.arange(float(args.gain_start),float(args.gain_stop)+float(args.gain_interval),float(args.gain_interval)))
gain_range = [ round(elem, 2) for elem in gain_range ]
res_range=list(numpy.arange(float(args.res_start),float(args.res_stop)+float(args.res_interval),float(args.res_interval)))
res_range = [ round(elem, 2) for elem in res_range ]
print "Resolution range is:"
print res_range
print "Gain range is:"
print gain_range
for g in gain_range:
	for r in res_range:
		name = args.cancer_name + '_' + args.metric + '_' + args.lens + '_' + str(r) + '_' + str(g)
		payload={"network_specifications": [ {"name": name,"column_set_id": args.column_set_id, "metric": {"id": metric} ,"lenses":[{"resolution":r,''"gain":g,"equalize":"false","id":len1},{"resolution":r,"gain":g,"equalize":"false","id":len2}]}]}
		'''payload={"network_specifications": [ {"name": name,"column_set_id": "-3616538532371804765", "metric": {"id": "Correlation"} ,"lenses":[{"resolution":r,''"gain":g,"equalize":"false","id":"Neighborhood Lens 1"},{"resolution":r,"gain":g,"equalize":"false","id":"Neighborhood Lens 2"}]}],"async":{}}'''
		print ('Creating network '+ str(int(g*r)) + ' of ' + str(len(gain_range)*len(res_range)) + ' : ' + name)
		a=session.post('https://core.ayasdi.com/v1/sources/'+args.source+'/networks',json.dumps(payload),headers=headers).content
		b=json.loads(a)
		net=b['id']
		dict[name]=net # Saving a dictionary of network names and ID's - Could be used later
		print('Downloading assigned network ID: ' + net)
		ParseAyasdiGraph(args.source,net,'uer2102@columbia.edu','ColumbiaAyasdi2015!',name)
'''print (dict)'''