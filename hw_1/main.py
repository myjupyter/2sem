#!/usr/bin/env python3

import download
import process

def main():
    with open('links.txt') as file:
        links = [link for link in file.read().split('\n') if len(link) != 0]

    paths = download.download_all(links)

    docs = process.Documents(paths)
    ss = docs.search_n('На чикагской выставке 1893 года посетители катались на колесе обозрения, а на выставке 1933 года — на скайрайде', 10)
    
    for s in ss:
        print('Текст: {:s}\nПредложение: {:s}\nВес: {:f}\n\n'.format(s.doc_name, s.__str__(), s.weight))

if __name__ == '__main__':
    main()
