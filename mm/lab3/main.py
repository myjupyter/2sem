#!/usr/bin/env python3

import parsers
import dictionary


def main():
    print(dictionary.bs_energy)
    mol = parsers.Molecule('mol.mol2')
    for k, b in mol.bonds.items():
        i, j = b.origin_atom_id, b.target_atom_id
        print(mol.atoms[i], mol.atoms[j])


if __name__ == '__main__':
    main()
