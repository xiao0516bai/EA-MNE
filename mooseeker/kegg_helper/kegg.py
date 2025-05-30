'''
Autor: caoyh
Date: 2022-11-20 11:33:15
LastEditTime: 2022-11-20 11:33:22
'''

import os
import re
import sys
import time
import json
import random
import requests
import threading
import queue
from time import sleep
from lxml import etree
from bs4 import BeautifulSoup
from functools import reduce
from rdkit import Chem
from concurrent.futures import ThreadPoolExecutor
import threading

# Define functions
def sWrite(string):
    sys.stdout.write(string)
    sys.stdout.flush()


def sError(string):
    sys.stderr.write(string)
    sys.stderr.flush()

def get_user_agent():
    file = 'KEGG/user_agent.txt'
    with open(file, 'r') as f:
        agents = [line.strip() for line in f]
    # for line in lines:
    #     agents.append(line.strip())
    return random.choice(agents)

def read_rxn_list():
    file = 'KEGG/rxn_list.txt'
    with open(file, 'r') as f:
        data = [line.strip() for line in f]
    return data

def read_cpd_list():
    file = 'KEGG/cpd_list.txt'
    with open(file, 'r') as f:
        data = [line.strip() for line in f]
    return data


## global variable
headers = {
    'User-Agent': get_user_agent()
}
# dynamic05.xingsudaili.com:10010

proxy = {
    # # # http://用户名:密码@代理的接口信息
    # "http": "http://bella:abc123@dynamic05.xingsudaili.com:10010",
    # "https": "http://bella:abc123@dynamic05.xingsudaili.com:10010"
}


def kegg_get(kegg_id, server="http://rest.kegg.jp"):
    """
    Downloads the raw REST text entry for the provided KEGG ID,
    via the KEGG rest API @ https://rest.kegg.jp/
    """
    # headers = {'User-Agent': get_user_agent()}
    # Construct query
    rest_address = "/".join([server,"get",kegg_id])
    # Set up con_attempts counter
    con_attempts = 0
    # Try connecting
    while con_attempts < 5:
        con_attempts += 1
        session = requests.Session()
        r = session.get(rest_address, headers=headers, proxies=proxy)
        # r = rget(rest_address, headers=headers)
        if r.status_code == 200:
            return r.text
        else:
            # Server not returning a result, try again
            time.sleep(2)
    # The connection attempt limit was reached
    sError("Warning: Unable to download KEGG data for '%s'.\n" % str(kegg_id))
    return None

def create_kegg_dict(kegg_text):
    """
    Parses a KEGG REST text record into a dictionary. Accepts a single record,
    as it stops after the first '///'. All lines under a key are presented as a
    list of lists.
    """
    # Initialize dictionary
    kegg_dict = {}
    # Iterate over lines in raw text
    for line in kegg_text.split("\n"):
        if line == "///":
            # The first entry ends here
            break
        line = re.split(" +", line.rstrip())
        if line[0] != "":
            key = line[0]
            line = line[1:]
            try:
                kegg_dict[key].append(line)
            except KeyError:
                kegg_dict[key] = [line]
        else:
            kegg_dict[key].append(line[1:])
    return kegg_dict

def kegg_smiles(kegg_id, server="http://rest.kegg.jp"):
    """Downloads a KEGG compound molecule object and converts it to SMILES."""
    # Ensure that the ID is a valid KEGG compound
    if not re.fullmatch("^C[0-9]{5}$", kegg_id):
        sError("\nWarning: '%s' is not a valid KEGG compound ID.\n" % str(kegg_id))
        return None
    # Set up the query
    rest_address = "/".join([server,"get","cpd:"+kegg_id,"mol"])
    # headers = {'User-Agent': get_user_agent()}
    # Set up con_attempts counter
    con_attempts = 0
    # Try connecting
    while con_attempts < 5:
        con_attempts += 1
        r = requests.get(rest_address, headers=headers, proxies=proxy)
        if r.status_code == 200:
            try:
                mol = Chem.MolFromMolBlock(r.text)
                smiles = Chem.MolToSmiles(mol)
                return smiles
            except:
                sError("\nWarning: SMILES could not be produced for KEGG ID '%s'.\n" % str(kegg_id))
                return None
        else:
            time.sleep(2)
    # The connection attempt limit was reached
    sError("\nWarning: Unable to download SMILES of molecule data for '%s'.\n" % str(kegg_id))
    return None


def threaded_kegg_get(queries):
    """Threaded implementation of kegg_get."""

    def worker():
        while True:
            query = work.get()
            if query is None:
                break
            # result = kegg_get(query)
            cpd = KEGGCompound(query)
            result = cpd.cpd_dict
            output.put((query, result))
            work.task_done()

    # Initialise queues
    work = queue.Queue()
    output = queue.Queue()

    for query in queries:
        work.put(query)

    # Start threads
    threads = []
    for i in range(16):
        t = threading.Thread(target=worker)
        t.start()
        threads.append(t)

    # Report on progress
    while True:
        if len(queries) == 0:
            progress = 100.0
        else:
            progress = float(output.qsize() / len(queries) * 100)
        sys.stdout.write("\rQuerying KEGG... %0.1f%%" % progress)
        sys.stdout.flush()
        if output.qsize() == len(queries):
            print("")
            break
        sleep(0.5)

    # Join work queue
    work.join()

    # Stop workers
    for i in range(len(threads)):
        work.put(None)
    for t in threads:
        t.join()

    # Get results
    results = {}
    while not output.empty():
        result = output.get()
        results[result[0]] = result[1]
    return results

class KEGGReaction(object):
    def __init__(self, kegg_id) -> None:
        self.kegg_dict = create_kegg_dict(kegg_get(kegg_id))

        eq_list = self.kegg_dict["EQUATION"][0]
        is_comp = re.compile("^C[0-9]{5}")

        self.rxn_dict = {
            'name': reduce(lambda x,y: x+y, self.kegg_dict['NAME']) if 'NAME' in self.kegg_dict.keys() else None,
            'kegg_id': kegg_id,
            'bigg_id': KeggId2BiggId(kegg_id),
            'equation': reduce(lambda x,y: x+' '+y, self.kegg_dict["EQUATION"][0]),
            'definition': reduce(lambda x,y: x+' '+y, self.kegg_dict["DEFINITION"][0]),
            'reactants': list(set(filter(is_comp.match, eq_list[0:eq_list.index("<=>")]))),
            'products': list(set(filter(is_comp.match, eq_list[eq_list.index("<=>")+1:]))),
            'enzyme': self.kegg_dict["ENZYME"][0] if "ENZYME" in self.kegg_dict.keys() else None
        }


class KEGGCompound(object):
    def __init__(self, kegg_id) -> None:
        self.kegg_dict = create_kegg_dict(kegg_get(kegg_id))

        self.cpd_dict = {
            'name': self._get_name(),
            'kegg_id': kegg_id,
            'bigg_id': KeggId2BiggId(kegg_id),
            'formula': self.kegg_dict["FORMULA"][0][0] if "FORMULA" in self.kegg_dict.keys() else None,
            'pubchem': self.getPubchemCid(),
            'chebi': self.getChebi(),
            'smile': kegg_smiles(kegg_id),
        }

    def _get_name(self):
        names = []
        for i in range(len(self.kegg_dict["NAME"])):
            if len(self.kegg_dict["NAME"][i]) > 1:          
                name = " ".join(self.kegg_dict["NAME"][i]).strip(";") 
            else:
                name = self.kegg_dict["NAME"][i][0].strip(";")
            names.append(name)      
        return names
            
    def getPubchemCid(self):        
        if "DBLINKS" in self.kegg_dict.keys(): 
            for ids in self.kegg_dict["DBLINKS"]:
                if ids[0] == 'PubChem:': return int(ids[1])  
        return None

    def getChebi(self):        
        if "DBLINKS" in self.kegg_dict.keys(): 
            for ids in self.kegg_dict["DBLINKS"]:
                if ids[0] == 'ChEBI:': return int(ids[1])
        return None
          
    def get_smile(self):
        return kegg_smiles(self.kegg_id)


def KeggId2BiggId(kegg_id):
    res = re.search('R|C', kegg_id).group()
    
    if res=='R':
        cookie = '_ga=GA1.2.208736496.1642561994; _gid=GA1.2.426276517.1656811686; _gat=1'
        database_source = 'kegg.reaction'

    elif res=='C':
        cookie = '_ga=GA1.2.208736496.1642561994; _gid=GA1.2.842153841.1656397739; _gat=1'
        database_source = 'kegg.compound'
    else:
        print('%s is a fault KEGG ID' %(kegg_id))
        return 
    
    url = 'http://bigg.ucsd.edu/advanced_search_external_id_results'

    headers = {
        'Cookie': cookie,
        'User-Agent': get_user_agent()}

    data = {
        'database_source': database_source,
        'query': kegg_id}

    con_attempts = 0
    while con_attempts < 5:
        con_attempts += 1
        r = requests.post(url=url, data=data, headers=headers, proxies=proxy)
        content = r.text
        if r.status_code == 200:
            details = etree.HTML(content)
            result = details.xpath('//div[3]/div/div/div/table/tbody/tr/td/a')
            if result != []: 
                bigg_ids = [res.text for res in result]
                return bigg_ids
            else:
                return None
        else:
            time.sleep(2)
    return None
    
def get_rxn_list():
    url = 'https://rest.kegg.jp/list/reaction/'
    r = requests.get(url=url, headers=headers, proxies=proxy)
    if r.status_code == 200:
        reactions = r.text
        reactions = reactions.replace('rn:', '').strip()
        reactions = reactions.split('\n')   
        rxn_list = [rxn.split('\t')[0] for rxn in reactions]
        return rxn_list
    else:
        print('Connection Error') 
        return None

def get_cpd_list():
    url = 'https://rest.kegg.jp/list/compound/'
    r = requests.get(url=url, headers=headers, proxies=proxy)
    if r.status_code == 200:
        compounds = r.text
        compounds = compounds.replace('cpd:', '').strip()
        compounds = compounds.split('\n')   
        cpd_list = [cpd.split('\t')[0] for cpd in compounds]
        return cpd_list
    else:
        print('Connection Error') 
        return None

def save_dict_json(file, data):
    with open(file, 'w') as f:
        json.dump(data, f)


