import re
import numpy

from scipy.optimize import minimize

from dictionary import bs_energy_K_r_mean, bs_energy_req_mean
from dictionary import bending_energy_K_q_mean, bending_energy_qeq_mean
from dictionary import atom_type_mapping, bs_energy, bending_energy

__token__ = '@<TRIPOS>'

def create_atoms(data=None):
    ss = data.split('\n')
    atoms = dict()
    for s in ss:
        if len(s) == 0:
            continue
        atom = Atom(s)
        atoms[atom.atom_id] = atom
    return atoms


def create_bonds(data=None):
    ss = data.split('\n')
    bonds = dict()
    for s in ss:
        if len(s) == 0:
            continue
        bond = Bond(s)
        bonds[bond.bond_id] = bond
    return bonds


__tokens_class_mapping__ = {
    'ATOM': create_atoms,
    'BOND': create_bonds,
}


class Molecule:
    def __init__(self, filepath):
        with open(filepath) as file:
            raw_text = file.read()
        sections = split_on_section(raw_text)
        obj = dict()
        for section_name, data in sections.items():
            if __tokens_class_mapping__.get(section_name) is not None:
                obj[section_name.lower()] = __tokens_class_mapping__[
                    section_name](data)
        self.__dict__.update(obj)
        for bond in self.bonds.values():
            bond.set_atoms_link(self.atoms)
        self.__bond_angle_count = -1

    def bond_stretching_energy(self):
        bse = 0
        for bond in self.bonds.values():
            i, j = bond.origin_atom_id, bond.target_atom_id
            rho = self.atoms[i].distance(self.atoms[j])
            x, name1, y, name2 = *self.atoms[i].atom_types, *self.atoms[j].atom_types
            K_r, req = bs_energy_K_r_mean, bs_energy_req_mean
            if bs_energy.get(frozenset((name1, name2))) is not None:
                coef = bs_energy.get(frozenset((name1, name2)))
                K_r, req = coef['K_r'], coef['req']
                coef = bs_energy[frozenset((name1, name2))]
            bse += K_r * (rho - req) ** 2 
        return bse / 2 

    def bse(self):
        rho_list = []
        for bond in self.bonds.values():
            i, j = bond.origin_atom_id, bond.target_atom_id
            rho_list.append(self.atoms[i].distance(self.atoms[j]))
        return self.__bond_stretching_energy(rho_list)

    def optimized_bse(self):
        n = len(self.bonds)
        return minimize(lambda x: numpy.sum(self.__bond_stretching_energy(x)), x0=[0]*n, method='SLSQP', bounds = [(0, None)]*n) 

    def optimized_be(self):
        self.__bond_triplets()
        m = self.__bond_angle_count
        return minimize(lambda x: numpy.sum(self.__bending_energy(x)), x0=[0]*m, method='SLSQP', bounds = [(1., 180.)]*m) 

    def optimized_energy(self):
        self.__bond_triplets()
        n = len(self.bonds)
        m = int(self.__bond_angle_count) 
        print(n, m)
        return minimize(
                lambda x: numpy.sum(
                    numpy.concatenate((
                        self.__bond_stretching_energy(x[:n]),
                        self.__bending_energy(x[n:])
                    ))
                ),
                x0=[0.01]*(n+m),
                method='BFGS',
                tol=1e-10,
                bounds=[(0., None)]*n + [(0., 180.)]*m
        )

    def __bond_stretching_energy(self, rho_list):
        bse = numpy.array([])
        for k, bond in enumerate(self.bonds.values()):
            i, j = bond.origin_atom_id, bond.target_atom_id
            name1, name2 = self.atoms[i].atom_mapped_type, self.atoms[j].atom_mapped_type
            K_r, req = bs_energy_K_r_mean, bs_energy_req_mean
            if bs_energy.get(frozenset((name1, name2))) is not None:
                coef = bs_energy.get(frozenset((name1, name2)))
                K_r, req = coef['K_r'], coef['req']
                coef = bs_energy[frozenset((name1, name2))]
            bse = numpy.append(bse, 0.5 * K_r * (rho_list[k] - req) ** 2) 
        return bse

    def bending_energy(self):
        be = 0
        for bond_i, bond_js in self.__bond_triplets().items():
            for bond_j in bond_js:
                q = angle(self.bonds[bond_i], self.bonds[bond_j]) 
                f, b = self.bonds[bond_i].chain(self.bonds[bond_j])
                key1 = frozenset([(i, self.atoms[j].atom_mapped_type) for i, j in enumerate(f)])
                K_q, qeq = bending_energy_K_q_mean, bending_energy_qeq_mean
                if bending_energy.get(key1) is None:
                    key2 = frozenset([(i, self.atoms[j].atom_mapped_type) for i, j in enumerate(b)])
                    if bending_energy.get(key2) is not None:
                        t = bending_energy[key2]
                        K_q, qeq = t['K_q'], t['qeq'] 
                else:
                    t = bending_energy[key1]
                    K_q, qeq = t['K_q'], t['qeq']  
                be += K_q * (q - qeq) ** 2
        return be / 2

    def __bending_energy(self, q):
        be = numpy.array([])
        index_i = 0
        for bond_i, bond_js in self.__bond_triplets().items():
            for bond_j in bond_js:
                f, b = self.bonds[bond_i].chain(self.bonds[bond_j])
                key1 = frozenset([(i, self.atoms[j].atom_mapped_type) for i, j in enumerate(f)])
                K_q, qeq = bending_energy_K_q_mean, bending_energy_qeq_mean
                if bending_energy.get(key1) is None:
                    key2 = frozenset([(i, self.atoms[j].atom_mapped_type) for i, j in enumerate(b)])
                    if bending_energy.get(key2) is not None:
                        t = bending_energy[key2]
                        K_q, qeq = t['K_q'], t['qeq'] 
                else:
                    t = bending_energy[key1]
                    K_q, qeq = t['K_q'], t['qeq']  
                be = numpy.append(be, 0.5 * K_q * (q[index_i] - qeq) ** 2)
                index_i += 1
        return be

    def __coloumb_energy(self):
        pass

    def bond_triplets(self):
        return self.__bond_triplets()

    def __bond_triplets(self):
        bonds_list = dict()
        for i, b1 in self.bonds.items():
            for j, b2 in self.bonds.items():
                if j <= i:
                    continue
                inter = b1.atom_indexes.intersection(
                        b2.atom_indexes)
                if inter != set():
                    if bonds_list.get(b1.bond_id) is None:
                       bonds_list[b1.bond_id] = list()
                    bonds_list[b1.bond_id].append(b2.bond_id)
        self.__bond_angle_count = 0
        for v in bonds_list.values():
            self.__bond_angle_count += len(v)
        return bonds_list
        
    @property
    def atoms(self):
        return self.atom

    @property
    def bonds(self):
        return self.bond


def split_on_section(raw_text):
    raw_text.split(__token__)
    sections = dict()

    for s in raw_text.split(__token__):
        i = s.find('\n')
        if i == -1:
            raise ValueError('wrong file')
        sections[s[:i].strip()] = s[i + 1:]
    return sections


class Atom:
    __data_order = {
        0: ('atom_id', int),
        1: ('atom_name', str),
        2: ('x', float),
        3: ('y', float),
        4: ('z', float),
        5: ('atom_type', str),
        6: ('subst_id', int),
        7: ('subst_name', str),
        8: ('charge', float),
    }

    def __init__(self, data=None):
        if data is None:
            raise ValueError('wrong data')
        ss = list(map(lambda x: x.strip(), filter(
            lambda x: len(x) != 0, re.split('\t| ', data))))
        if len(ss) < len(Atom.__data_order):
            raise ValueError('not enough parameters')

        for i, v in enumerate(ss[:len(Atom.__data_order)]):
            name, func = Atom.__data_order[i]
            self.__dict__[name] = func(v)

    def __repr__(self):
        return 'Atom(id={}, name={}, charge={})'.format(
            self.atom_id, self.atom_name, self.charge)

    def __str__(self):
        return 'id={}, name={}, charge={}'.format(
            self.atom_id, self.atom_name, self.charge)

    def distance(self, atom = None) -> float:
        if atom is None:
            raise ValueError('None atom')
        return numpy.linalg.norm(self.coord - atom.coord)

    @property
    def atom_types(self):
        return (self.atom_type, self.atom_mapped_type)

    @property
    def atom_mapped_type(self):
        return atom_type_mapping.get(self.atom_type)
        
    @property
    def coord(self):
        return numpy.array([
            self.x,
            self.y,
            self.z,
        ])


class Bond:
    __data_order = {
        0: ('bond_id', int),
        1: ('origin_atom_id', int),
        2: ('target_atom_id', int),
        3: ('bond_type', str),
    }

    def __init__(self, data=None):
        if data is None:
            raise ValueError('wrong data')
        ss = list(map(lambda x: x.strip(), filter(
            lambda x: len(x) != 0, re.split('\t| ', data))))
        if len(ss) < len(Bond.__data_order):
            raise ValueError('not enough parameters')

        for i, v in enumerate(ss[:len(Bond.__data_order)]):
            name, func = Bond.__data_order[i]
            self.__dict__[name] = func(v)

    @property
    def atom_indexes(self):
        return set((self.origin_atom_id, self.target_atom_id))

    def get_atom_neighbour(self, atom_id):
        ids = self.atom_indexes
        if atom_id not in ids:
            return None
        return list(self.atom_indexes.difference(set([atom_id])))[0]

    def chain(self, bond):
        inter = self.intersect(bond)
        if inter == set():
            return None
        l = list(inter)[0]
        p = self.get_atom_neighbour(l)
        q = bond.get_atom_neighbour(l)
        t = (p, l, q)
        return tuple([[e for e in t[::j]] for j in [1,-1]])  

    def intersect(self, bond):
        return self.atom_indexes.intersection(bond.atom_indexes) 

    def set_atoms_link(self, atoms):
        self.atoms_link = atoms

    def vectorize(self):
        return self.atoms_link[self.origin_atom_id].coord - self.atoms_link[self.target_atom_id].coord

def cosine(x, y):
    return numpy.dot(x, y) / (numpy.linalg.norm(x) * numpy.linalg.norm(y)) 

def angle(x, y):
    d = numpy.degrees(numpy.arccos(cosine(x.vectorize(), y.vectorize())))
    return d if d > 90 else 180 - d

