#!/usr/bin/env python3

import sys
import os
import argparse
import json
import math

class Args:
    pass

def dcg(rel_list):
    return sum([(r / math.log2(i + 2)) for i, r in enumerate(rel_list)])

class NDCG:
    def __init__(self, docs, result):
        with open(result) as file:
            self.__rsentences = json.load(fp=file)['search_result']
        
        self.__docs = {}
        for doc in docs:
            with open(doc) as file:
                ls = json.load(fp=file)['data']
                for l in ls:
                    self.__docs[l['sentence']] = l['rel']


    def idcg(self, n = 1):
        return dcg(list(map(lambda x: x[1], sorted(self.__docs.items(), key=lambda x: x[1], reverse=True)))[:n]) 

    def ndcg(self):
        rel_val = list(map(lambda x: 0 if x is None else x, [self.__docs.get(s['sentence']) for s in self.__rsentences]))
        return dcg(rel_val) / self.idcg(len(rel_val))

def main():
    parser = argparse.ArgumentParser(description='nDCG') 
    parser.add_argument(
        '-d', '--docs',
        nargs = '+',
        action = 'append',
        required = True,
        help = 'template documents'
    )
    parser.add_argument(
        '-r', '--results',
        type = str,
        required = True,
        help = 'result document'
    )
    args = parser.parse_args(args=sys.argv[1:], namespace=Args())
    
    n = NDCG([''.join(d) for d in args.docs], args.results)
    print(n.ndcg())

if __name__ == '__main__':
    main()
