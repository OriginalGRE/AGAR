"""
Exports the LammpsDataWriter class
Could probably be implemented more elegantly as a subclass of some generic data writer; will have to figure this out
"""

from Lattice import Lattice


class LammpsDataWriter:

    DEFAULT_FILENAME = 'AGAR_out.data'

    def __init__(self, lattice: Lattice, filename: str = DEFAULT_FILENAME):
        """
        Accepts a Lattice object and outputs a corresponding LAMMPS data file
        :param lattice: Name of Lattice object
        :param filename: String; filename to save data file to
        """
        self.atoms = lattice.full_atom_list
        self.bonds = lattice.full_bond_list
        self.filename = filename

        boxmin, boxmax = [], []
        for i in range(3):
            boxmax.append(max(atom.position[i]+1.0 for atom in self.atoms))
            boxmin.append(min(atom.position[i]-1.0 for atom in self.atoms))

        header = ['LAMMPS data file through AGAR\n',
                  '\n',
                  '{} atoms\n'.format(len(self.atoms)),
                  '2 atom types\n',
                  '{} bonds\n'.format(len(self.bonds)),
                  '1 bond types\n']

        dims = ['{} {} xlo xhi\n'.format(boxmin[0], boxmax[0]),
                '{} {} ylo yhi\n'.format(boxmin[1], boxmax[1]),
                '{} {} zlo zhi\n'.format(boxmin[2], boxmax[2])]

        masses = ['Masses\n',
                  '\n',
                  '1 1\n',
                  '2 1\n']

        atom_section = ['Atoms # full\n',
                        '\n']

        for atom in self.atoms:
            if atom.owning_cell:
                atom_type = 1
            else:
                atom_type = 2
            atom_section.append('{} 1 {} 0 {} {} {} 0 0 0\n'.format(atom.id, atom_type, atom.position[0], atom.position[1], atom.position[2]))

        bond_section = ['Bonds\n',
                        '\n']
        for bond in self.bonds:
            bond_section.append('{} 1 {} {}\n'.format(bond['bond_id'], bond['point1'], bond['point2']))

        data_file_lines = header + ['\n'] + dims + ['\n'] + masses + ['\n'] + atom_section + ['\n'] + bond_section

        with open(self.filename, 'w') as datafile:
            datafile.writelines(data_file_lines)
