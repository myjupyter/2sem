import os
import requests
import re
import string

from urllib.parse import urlparse, unquote
from bs4 import BeautifulSoup

def __mkdirp(path_to):
    try:
        os.makedirs(path_to)
    except FileExistsError:
        pass

def __clear(text):
    soup = BeautifulSoup(text, features='html5lib')
    div = soup.findAll('div', {'id':'bodyContent'})
    return ' '.join([p.get_text().strip() for p in soup.select('p')])

def download_all(links):
    files = [] 
    for link in links:
        files.append(download(link))
    
    return files

def download(link):
    url = urlparse(unquote(link))
    
    file_name = os.path.basename(url.path)
    path = url.netloc+os.path.dirname(url.path)
    
    full_path = os.path.join(path, file_name)

    __mkdirp(path) 

    res = requests.get(link)

    if res.status_code != 200:
        raise Exception('{:s} has {:d} status code'.format(link, res.status_code))
    
    text = __clear(res.text)
    with open(full_path, 'w') as file:
        file.write(text)
   
    return full_path
