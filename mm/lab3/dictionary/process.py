def process_data(data : dict):
    new_data = dict()
    new_list = list()
    for el in data["data"]:
        el["bind"] = ''.join(el["bind"].split(' '))
        new_list.append(el)
    return {"data": new_list}
