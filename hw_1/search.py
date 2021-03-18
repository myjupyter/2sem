#!/usr/bin/env python3

import sys
import json

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
    
    rr = []
    with open(args.output, 'w') as file:
        for s in ss:
            r = {
                'text': s.doc_name,
                'sentence': s.__str__(),
                'weight': s.weight,
            }
            rr.append(r)
        o = {
            'query': args.query,
            'links': links,
            'search_result': rr,
        }
        json.dump(o, file, ensure_ascii=False,indent=4)

if __name__ == '__main__':
    main()
