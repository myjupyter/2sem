import re
import numpy

__bs_energy_bond_sep_regex = re.compile(' |-')

def process_data(data : dict):
    new_data = dict()
    new_list = list()
    for el in data["data"]:
        el["bind"] = ''.join(el["bind"].split(' '))
        new_list.append(el)
    return {"data": new_list}

def process_bs_energy(data: dict):
    new_data = dict()
    for el in data["data"]:
        key = frozenset(filter(None, __bs_energy_bond_sep_regex.split(el['bind'])))
        del el['bind']
        new_data[key] = el 
    return new_data

def K_r_mean(data: dict):
    K_r_list = list()
    for el in data['data']:
        K_r_list.append(el['K_r'])
    return numpy.mean(K_r_list)

def req_mean(data: dict):
    req_list = list()
    for el in data['data']:
        req_list.append(el['req'])
    return numpy.mean(req_list)
