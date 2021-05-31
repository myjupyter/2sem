import re
import numpy

from itertools import combinations
from scipy.optimize import minimize

from dictionary import bs_energy_K_r_mean, bs_energy_req_mean
from dictionary import bending_energy_K_q_mean, bending_energy_qeq_mean
from dictionary import atom_type_mapping, bs_energy, bending_energy, van_der_waals_energy

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

class Indexer:
    def __init__(self, xk, n = None):
        if n is None or len(xk) - 3 != 3 * n:
            raise ValueError('wrong atoms number: got {} want {}'.format(len(xk), 3 * n))
        self.xk = xk
    
    def __getitem__(self, i):
        return numpy.array(self.xk[3*i:3*i+3])
 
class Parameters:
    def __init__(self, n = None, bonds = None, triplets = None):
        self.n = n
        self.bonds = bonds
        self.triplets = triplets

def bond_stretching_energy(parameters):
    def energy(xk):
        res = numpy.array([])
        ind = Indexer(xk, parameters.n)
        for bond in parameters.bonds:
            r = numpy.linalg.norm(ind[bond.origin_atom_id] - ind[bond.target_atom_id])
            res = numpy.append(res, 0.5 * bond.K_r * (bond.req - r) ** 2)
        return res
    return energy 


def bending_energy(parameters):
    def energy(xk):
        res = numpy.array([])
        ind = Indexer(xk, parameters.n)
        for triplet in parameters.triplets:
            #q = ...
            res = numpy.append(res, 0.5 * triplet.K_q * (triplet.qeq - q) ** 2) 
        return res
    return energy

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

    def get_ordered_coords():
        return [atom.coord for _, atom in sorted(self.atoms.items(), key=lambda x: x[0])]

    def bond_stretching_energy(self):
        coords = self.get_ordered_coords()
        xk = numpy.array([])
        for coord in coords:
            xk = numpy.append(xk, coord)
        return numpy.sum(self.__bond_stretching_energy(xk))

    def __bond_stretching_energy(self, xk):
        nulls = numpy.array([0.]*3)
        return bond_stretching_energy(
                Parameters(
                    n = len(self.atoms),
                    bonds = self.bonds.values(),
                )
        )(numpy.append(nulls, xk))

    def optimized_energy(self):
        def func(xk):
            return numpy.sum(
                self.__bond_stretching_energy(xk),
            )
        n = len(self.atoms) * 3
        return minimize(
            func,
            x0=numpy.fabs(numpy.random.rand(n)) * 100,
            method='SLSQP',
            options = {
                'maxiter': 2000,
                'ftol': 10e-10,
            },
            bounds=[(0., None)]*n, 
            constraints = [
                {
                    'type': 'eq', 
                    'fun': func,
                }
            ]
        )
 
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

    def end_combinations(self):
        bonds_list = self.__bond_triplets()
        ec = set()
        for b1, blist in bonds_list.items():
            for b2 in blist:
                e = self.bonds[b1].chain(self.bonds[b2])[0]
                ec.add((e[0], e[2]))

        for bond in self.bonds.values():
            ec.add(bond.atom_indexes_tuple)

        return ec
    
    def van_der_waals_atom_pairs(self):
        all_pairs = set(combinations(self.atoms.keys(), r=2))
        return all_pairs.difference(self.end_combinations())

    def optimized_van_der_waals_energy(self):
        p = len(self.van_der_waals_atom_pairs())
        return minimize(
                lambda x: numpy.sum(
                    numpy.concatenate((
                        self.__van_der_waals_energy(x),
                    ))
                ),
                x0=[10.0]*p,#numpy.fabs(numpy.random.rand(p)) * 100,
                method='SLSQP',
                options= {
                    'maxiter': 100000,
                    'ftol': 10e-20,
                },
                bounds=[(0., None)]*p,
                constraints=[{'type': 'eq', 'fun': lambda x: numpy.sum(self.__van_der_waals_energy(x))}],
        )



    def __van_der_waals_energy(self, rho_list):
        seq = numpy.array([]) 
        for i, t in enumerate(self.van_der_waals_atom_pairs()):
            f, s = t
            atom1, atom2 = self.atoms[f], self.atoms[s]
            R_ij = atom1.R + atom2.R
            e_ij = (atom1.e() * atom2.e()) ** 0.5
            A = e_ij * (R_ij ** 12)
            B = 2 * e_ij * (R_ij ** 6)
            seq = numpy.append(seq,  A / rho_list[i] ** 12 - B / rho_list[i] ** 6) 
        return seq
    
    def __van_der_waals_and_coloumb_energy(self, rho_list):
        seq = numpy.array([]) 
        eps = 8.9875517873681764*(10**9)
        for i, t in enumerate(self.van_der_waals_atom_pairs()):
            f, s = t
            atom1, atom2 = self.atoms[f], self.atoms[s]
            R_ij = atom1.R + atom2.R
            e_ij = (atom1.e * atom2.e) ** 0.5
            A = e_ij * (R_ij ** 12)
            B = 2 * e_ij * (R_ij ** 6)
            seq = numpy.append(seq,  A / (rho_list[i] ** 12) - B / (rho_list[i] ** 6) + eps * atom1.charge * atom2.charge * 10**(-18) / rho_list[i]) 
        return seq

    def van_der_waals_energy(self):
        ss = 0
        for f, s in self.van_der_waals_atom_pairs():
            atom1, atom2 = self.atoms[f], self.atoms[s]
            R = atom1.distance(atom2)
            R_ij = atom1.R + atom2.R
            e_ij = (atom1.e * atom2.e) ** 0.5
            A = e_ij * (R_ij ** 12)
            B = 2 * e_ij * (R_ij ** 6)
            ss += A / (R ** 12) - B / (R ** 6) 
        return ss 

    def coloumb_energy(self):
        ss = 0
        eps = 8.9875517873681764*(10**9)
        for f, s in self.van_der_waals_atom_pairs():
            atom1, atom2 = self.atoms[f], self.atoms[s]
            R = atom1.distance(atom2)
            ss += eps * atom1.charge * atom2.charge * 10**(-18) / R
        return ss
        
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
    def R(self):
        return self.__vdw_coef('R*j')

    @property
    def e(self):
        return self.__vdw_coef('e_k') 

    def __vdw_coef(self, coef_name):
        return van_der_waals_energy.get(self.atom_mapped_type).get(coef_name)

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

class Triplets:
    def __init__(self, b1, b2):
        self.b1 = b1
        self.b2 = b2

    def __coef(self):
        // 
        return coef if coef is not None else {'K_q': bending_energy_K_q_mean, 'qeq': bending_energy_qeq_mean} 

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
            
    def __repr__(self):
        return 'Bond(id={}, origin_atom_id={}, target_atom_id={})'.format(
            self.bond_id, self.origin_atom_id, self.target_atom_id)

    def __coef(self):
        x, y = self.atoms_link[self.origin_atom_id].atom_mapped_type, self.atoms_link[self.target_atom_id].atom_mapped_type 
        coef = bs_energy.get(frozenset((x, y)))
        return coef if coef is not None else {'K_r': bs_energy_K_r_mean, 'req': bs_energy_req_mean} 

    @property
    def K_r(self):
        return self.__coef()['K_r']

    @property
    def req(self):
        return self.__coef()['req']

    @property
    def atom_indexes_tuple(self):
        return (self.origin_atom_id, self.target_atom_id)

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

def anglee(x, y):
    d = numpy.degrees(numpy.arccos(cosine(x,y)))
    return d if d > 90 else 180 - d
