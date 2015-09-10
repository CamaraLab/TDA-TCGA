import json
import requests

headers = {"Content-type": "application/json"}
session = requests.Session()
session.post('https://core.ayasdi.com/login', data={'username': 'uer2102@columbia.edu', 'passphrase': 'ColumbiaAyasdi2015!'})

res=[51,52,53,54]
gain=[1,2,3,4]
for g in gain:
	for r in res:
		payload={"network_specifications": [ {"name": str(r)+'_'+str(g),"column_set_id": "-3616538532371804765", "metric": {"id": "Correlation"} ,"lenses":[{"resolution":r,"gain":g,"equalize":"false","id":"Neighborhood Lens 1"},{"resolution":r,"gain":g,"equalize":"false","id":"Neighborhood Lens 2"}]}]}
		session.post('https://core.ayasdi.com/v1/sources/1440532486038/networks',json.dumps(payload),headers=headers)



"""
payload={"network_specifications": [ {"name": "yolo","column_set_id": "-3616538532371804765", "metric": {"id": "Correlation"} ,"lenses": [ {"id": "Mean","resolution": 20, "gain": 2.5, "equalize": "TRUE"} ]}]}
session.post('https://core.ayasdi.com/v1/sources/1440532486038/networks',json.dumps(payload),headers=headers).content


payload={"network_specifications": [ {"name": "yolo1","column_set_id": "-3616538532371804765", "metric": {"id": "Correlation"}, "metric": {"id": "Correlation"} ,"lenses": [ {"id": "Mean","resolution": 20, "gain": 2.5, "equalize": "TRUE"} ],"lenses": [ {"id": "Mean","resolution": 20, "gain": 2.5, "equalize": "TRUE"} ]}]}
session.post('https://core.ayasdi.com/v1/sources/1440532486038/networks',json.dumps(payload),headers=headers).content

"lenses":[{"resolution":36,"gain":2.999999999999999,"equalize":false,"id":"Neighborhood Lens 1"},{"resolution":37,"gain":2.999999999999999,"equalize":false,"id":"Neighborhood Lens 2"}]

{
   network_specifications: [ {
     name: string,
     column_set_id: string,
     metric: {
       id: string,
       source_id: string
     } ,
     lenses: [ {
       id: string,
       column_index: long,
       resolution: long,
       gain: double,
       equalize: Boolean
     } ]
   } ],
   name: string,
   column_set_id: string,
   metric: {
     id: string,
     source_id: string
   } ,
   lenses: [ {
     id: string,
     column_index: long,
     resolution: long,
     gain: double,
     equalize: Boolean
   } ],
   async: {
     callback: URI
   } 
}

"""