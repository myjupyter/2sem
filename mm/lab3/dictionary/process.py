import re
import numpy

__bind_sep_regex = re.compile(' |-')

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
        key = frozenset(filter(None, __bind_sep_regex.split(el['bind'])))
        del el['bind']
        new_data[key] = el 
    return new_data

def process_bending_energy(data: dict):
    new_data = dict()
    for el in data['data']:
        atoms = list(filter(None, __bind_sep_regex.split(el['bind']))) 
        del el['bind']
        new_data[frozenset([(i, atom) for i, atom in enumerate(atoms)])] = el
        new_data[frozenset([(i, atom) for i, atom in enumerate(atoms[::-1])])] = el
    return new_data

def K_q_mean(data: dict):
    K_q_list = list()
    for el in data['data']:
        K_q_list.append(el['K_q'])
    return numpy.mean(K_q_list)

def K_r_mean(data: dict):
    K_r_list = list()
    for el in data['data']:
        K_r_list.append(el['K_r'])
    return numpy.mean(K_r_list)

def qeq_mean(data: dict):
    qeq_list = list()
    for el in data ['data']:
        qeq_list.append(el['qeq'])
    return numpy.mean(qeq_list)

def req_mean(data: dict):
    req_list = list()
    for el in data['data']:
        req_list.append(el['req'])
    return numpy.mean(req_list)
