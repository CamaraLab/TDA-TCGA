
# coding: utf-8

# In[191]:

import networkx
import simplejson as json
import numpy
import numpy.linalg
import random
import requests
import pickle
import csv
import os
import pandas as pd


# In[192]:

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


# In[216]:

import argparse
parser = argparse.ArgumentParser()
parser.add_argument("cancer_name")
parser.add_argument("file1")
args = parser.parse_args()

#file1='C:\Users\Udi\SkyDrive\TCGA_CURATED\OTHER_CUR\KIRC_CUR\KIRC_CUR_BIG_matrix4.csv'
file1=args.file1
print "Importing matrix file :" + file1
df=pd.read_csv(file1,index_col=0)


# In[176]:

exp=df.filter(like="exp_")
b=exp.var(axis=0)
b=b.sort_values(ascending=0)
top4500=list(b[:4500].index)
print "Done"
print b[:4500]

# In[169]:

print "Establishing connection to aysdi"
from ayasdi.core.api import Api
connection = Api(username="uer2102@columbia.edu", password="ColumbiaAyasdi2015!")


# In[170]:

print "Uploading matrix file"
source=connection.upload_source(file1)


# In[177]:

print "Generating column set of top 4500 genes"
col_set = source.create_column_set(top4500, "exp_top4500")

print "Source name: "+ source.id
print "Column set: " + col_set['id']
#l1=pd.DataFrame(col_set)
#l1


# In[180]:

#Skelaton network
a={u'column_set_id': col_set['id'],
 u'lenses': [{u'equalize': False,
   u'gain': 3.0000000000000004,
   u'id': u'Neighborhood Lens 1',
   u'resolution': 30},
  {u'equalize': False,
   u'gain': 3.0000000000000004,
   u'id': u'Neighborhood Lens 2',
   u'resolution': 30}],
 u'metric': {u'id': u'Correlation'}}


# In[182]:

#network = source.create_network("test40", a)


# In[195]:

#gain_range=list(numpy.arange(float(args.gain_start),float(args.gain_stop)+float(args.gain_interval),float(args.gain_interval)))
#gain_range = [ round(elem, 2) for elem in gain_range ]
#res_range=list(numpy.arange(float(args.res_start),float(args.res_stop)+float(args.res_interval),float(args.res_interval)))
#res_range = [ int(elem) for elem in res_range ]

cancer_name=args.cancer_name # arg.cancer
res_range=[10,20,30,40,50,60,70,80] # if range
gain_range=[1.5,2.5,3.5,4.5,5.5,6.5,7.5,8.5] # if gain 
counter=0


# In[212]:

for g in gain_range:
    for r in res_range:
        counter+=1
        name = cancer_name + '_' + "Cor" + '_' + "Neigh" + '_' + str(r) + '_' + str(g)
        if os.path.isfile(name+'.json') and os.path.isfile(name+'.gexf'):
            print (name + " network exists, skipping")
        else:
			a['lenses'][0]['gain']=g
			a['lenses'][1]['gain']=g
			a['lenses'][0]['resolution']=r
			a['lenses'][1]['resolution']=r
			print ('Creating network '+ str(counter) + ' of ' + str(len(gain_range)*len(res_range)) + ' : ' + name)
			net = source.create_network(name, a)
            #print net['id']
            #if args.download
			ParseAyasdiGraph(source.id,net.id,'uer2102@columbia.edu','ColumbiaAyasdi2015!',name)