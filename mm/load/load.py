import os
import requests

def mkdirp(path: str):
    os.makedirs(path, exist_ok = True)

class Loader:
    def __init__(self):
        pass

    def download_to(self, url: str, path):
        mkdirp(path)

        r = requests.get(url)
        if r.status_code != 200:
            raise Exception('Status {}'.format(r.status_code))

        file_name = os.path.basename(url)
        with open(os.path.join(path,file_name), 'w') as file:
            file.write(r.text)

    r 
