#!/usr/bin/env python3

import sys
import os
import argparse
import json

class Args:
    pass

class Template:
    def __init__(self, filepath = None):
        with open(filepath) as file:
            self.__sentences = file.read().split('.') 
        self.__filename = os.path.basename(filepath)

    def write(self, path):
        data = {'data' : []}
        for sentence in self.__sentences:
            data['data'].append({'rel' : 0, 'sentence': sentence})

        with open(os.path.join(path, self.__filename + '.json'), 'w') as file:
            json.dump(data, fp=file, indent=4, ensure_ascii=False)

def main():
    parser = argparse.ArgumentParser(description='template maker')
    parser.add_argument(
        '-p', '--path',
        required = True,
        type = str,
        help = 'path with raw files',
    )
    parser.add_argument(
        '-o', '--output',
        type = str,
        default = './',
        help = 'path for new templates',
    )
    args = parser.parse_args(args=sys.argv[1:], namespace=Args())

    files_list = [os.path.join(args.path, f) for f in os.listdir(args.path) if os.path.isfile(os.path.join(args.path, f))]
    for filepath in files_list:
        t = Template(filepath)
        t.write(args.output)

if __name__ == '__main__':
    main()
