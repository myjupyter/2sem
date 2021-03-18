#!/usr/bin/env python3

import sys

import download
import process
import flags

def main():
    parser = flags.ArgsParser()
    args = parser.Parse(sys.argv[1:])

    with open(args.links) as file:
        links = [link for link in file.read().split('\n') if len(link) != 0]

    paths = download.download_all(links, args.download)

    docs = process.Documents(paths)
    ss = docs.search_n(args.query, args.number, args.method)
    
    for s in ss:
        print('Текст: {:s}\nПредложение: {:s}\nВес: {:f}\n\n'.format(s.doc_name, s.__str__(), s.weight))

if __name__ == '__main__':
    main()
