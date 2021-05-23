#!/usr/bin/env python3

import parsers
import dictionary


def main():
    mol = parsers.Molecule('mol.mol2')
    print(mol.bond_stretching_energy())
    print(mol.bending_energy())

if __name__ == '__main__':
    main()
