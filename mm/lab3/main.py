#!/usr/bin/env python3

import parsers
import dictionary


def main():
    mol = parsers.Molecule('mol.mol2')
    # for k, b in mol.bonds.items():
    #     i, j = b.origin_atom_id, b.target_atom_id
    #     print(mol.atoms[i], mol.atoms[j], f'distance is {mol.atoms[i].distance(mol.atoms[j])}')
    print(mol.bond_stretching_energy())

if __name__ == '__main__':
    main()
