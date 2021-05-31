#!/usr/bin/env python3

import parsers
import dictionary


def main():
    mol = parsers.Molecule('test.mol2')
    #print(mol.bond_stretching_energy())
    # print(mol.van_der_waals_energy())
    #print(mol.van_der_waals_atom_pairs())
    #print(mol.bse())
    # print(mol.optimized_be())
    # print(mol.optimized_van_der_waals_energy())
    # print(mol.optimized_energy())
    # print(mol.coloumb_energy())
    #print(mol.bond_stretching_energy())
    print(mol.optimized_energy())

if __name__ == '__main__':
    main()
