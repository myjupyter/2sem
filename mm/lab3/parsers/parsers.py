import re
import numpy

__token__ = '@<TRIPOS>'

def create_atoms(data = None):
    ss = data.split('\n')
    atoms = dict()
    for s in ss:
        if len(s) == 0:
            continue
        atom = Atom(s)
        atoms[atom.atom_id] = atom
    return atoms

def create_bonds(data = None):
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
                obj[section_name.lower()] = __tokens_class_mapping__[section_name](data)        
        self.__dict__.update(obj)
    
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
        sections[s[:i].strip()] = s[i+1:]
    return sections




class Atom:
    __data_order = {
            0: 'atom_id',
            1: 'atom_name',
            2: 'x',
            3: 'y',
            4: 'z',
            5: 'atom_type',
            6: 'subst_id',
            7: 'subst_name',
            8: 'charge',
    }

    def __init__(self, data = None):
        if data is None:
            raise ValueError('wrong data')
        ss = list(map(lambda x: x.strip(), filter(lambda x: len(x) != 0, re.split('\t| ', data) )))
        if len(ss) < len(Atom.__data_order):
            raise ValueError('not enough parameters')

        for i, v in enumerate(ss[:len(Atom.__data_order)]):
            self.__dict__[Atom.__data_order[i]] = v

    def __repr__(self):
        return 'Atom(id={}, name={}, charge={})'.format(self.atom_id, self.atom_name, self.charge)

    def __str__(self):
        return 'id={}, name={}, charge={}'.format(self.atom_id, self.atom_name, self.charge)

    @property
    def coord(self):
        return numpy.array([
            self.x,
            self.y,
            self.z,
        ])
    

class Bond:
    __data_order = {
        0: 'bond_id',
        1: 'origin_atom_id',
        2: 'target_atom_id',
        3: 'bond_type',
    }

    def __init__(self, data = None):
        if data is None:
            raise ValueError('wrong data')
        ss = list(map(lambda x: x.strip(), filter(lambda x: len(x) != 0, re.split('\t| ', data))))
        if len(ss) < len(Bond.__data_order):
            raise ValueError('not enough parameters')

        for i, v in enumerate(ss[:len(Bond.__data_order)]):
            self.__dict__[Bond.__data_order[i]] = v
