#!/usr/bin/env python3

import rdkit

import load

url = 'http://130.92.106.217:8080/GDB-13_Subsets/GDB13_Subset-ABCDEFGH.smi.gz'

def main():
    l =  load.Loader()
    l.download_to(url, './default/GDB13')


if __name__ == '__main__':
    main()
