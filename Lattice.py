"""
Exports the Lattice class for AGAR
"""
from typing import Optional, Tuple
from Cell import Cell, Connection
from itertools import product
import numpy as np


class Atom:
    """
    Data structure for use in Lattice.
    Manual addition of atoms to the system is not recommended; if at all possible, edit the cell structure instead.
    """
    __slots__ = ["id", "position", "point_type", "owning_cell", "bonds"]

    def __init__(self, atom_id: int, position: Tuple[float, float, float], point_type: Optional[int] = None,
                 owning_cell: Optional[Tuple[int, int, int]] = None):
        """
        Defines the atom for atom list creation
        :param atom_id: Atom identifier; uniqueness is assumed but not not enforced
        :param position: Atom position in (x, y, z) format
        :param point_type: Optional; point in Cell structure the atom is associated with. Used to set up connectivity.
        :param owning_cell: Optional; specific unit cell the atom is associated with. Used to set up connectivity.
        """
        self.id = atom_id
        self.position = position
        self.point_type = point_type
        self.owning_cell = owning_cell
        self.bonds = 0

class Lattice:
    """
    Provides methods to generate atom and bond lists for a full lattice from a Cell object
    Cells are assumed to originate at (0, 0, 0) and extend in positive direction
    Deviations from this will probably work, but may result in offsets in atom positions
    """

    DEFAULT_SCALE = 1
    ALLOW_DANGLING = False
    MIRROR_DANGLING = True
    KILL_DANGLING = False
    def __init__(self, cell: Cell, size: tuple[int, int, int], scale: Optional[float] = DEFAULT_SCALE,
                 dangling: Optional[bool] = ALLOW_DANGLING, symmetric: Optional[bool] = MIRROR_DANGLING, prune: Optional[bool] = KILL_DANGLING ):
        """
        Constructor; calls self.build_lists to create atom and bond lists from Cell object
        :param cell: Cell object to be tiled
        :param size: Lattice size in terms of unit cells per dimension; e.g. (1, 1, 1) represents a single cell
        :param scale: Factor to scale cell dimensions and point positions by
        :param dangling: Whether bonds outside the current cell are constructed
        :param symmetric: Whether dangling bonds are made symmetric across the lattice via mirrored bonds at edge cells
        :param prune: Whether points with only one connection are recursively pruned out of the connection list
        Note that (dangling = True, prune = True) is not equivalent to (dangling = False)
        """
        self.lattice_size = size

        self.base_point_list = cell.points
        for point in self.base_point_list:
            point.coords = np.multiply(point.coords, scale)

        self.base_connection_list = cell.connections
        self.cell_dims = np.multiply(cell.dims, scale)
        self.dangling = dangling
        self.symmetric = symmetric
        self.prune = prune

        size_x, size_y, size_z = range(self.lattice_size[0]), range(self.lattice_size[1]), range(self.lattice_size[2])
        self.cell_list = list(product(size_x, size_y, size_z))

        self.vector_list, self.full_atom_list, self.full_bond_list = [], [], []
        self.atom_count, self.bond_count = 0, 0
        self.float_safety = 1E-6
        self.build_lists()

    def build_lists(self):
        """
        Generates atom and bond lists for the current Lattice instance
        Separated from the constructor to enable more advanced editing; e.g. irregular cell lists

        First pass over all cells generates all points and (if enabled) a list of dangling bonds
        Second pass looks up all bonds attached to the current cell and populates them with atoms
        """
        # Dislike parts of this implementation; still looking to improve

        size_x, size_y, size_z = range(self.lattice_size[0]), range(self.lattice_size[1]), range(self.lattice_size[2])
        cell_list = list(product(size_x, size_y, size_z))
        mirrored_bond_list = []
        for current_cell in cell_list:
            current_cell_origin = np.multiply(current_cell, self.cell_dims)

            for point in self.base_point_list:
                self.atom_count += 1
                self.full_atom_list.append(Atom(self.atom_count, tuple(np.add(current_cell_origin, point.coords)),
                                                point.point_id, current_cell))

            if self.dangling is True:
                for connection in self.base_connection_list:
                    target_cell = tuple(np.add(current_cell, connection.cell_crossings))
                    target_cell_origin = np.multiply(target_cell, self.cell_dims)

                    if target_cell not in cell_list:
                        partner_coords = next(item.coords for item in self.base_point_list
                                              if item.point_id == connection.partner)
                        if not next((atom_check for atom_check in self.full_atom_list if (atom_check.point_type == connection.partner and atom_check.owning_cell == target_cell)), None):
                            self.atom_count += 1
                            self.full_atom_list.append(Atom(self.atom_count, tuple(np.add(target_cell_origin, partner_coords)),
                                                        connection.partner, target_cell))

                    if self.symmetric is True:
                        mirrored_target_cell = tuple(np.subtract(current_cell, connection.cell_crossings))
                        if mirrored_target_cell not in cell_list:
                            mirrored_cell_origin = np.multiply(mirrored_target_cell, self.cell_dims)
                            pseudo_origin_coords = next(
                                item.coords for item in self.base_point_list if item.point_id == connection.origin)
                            if not next((atom_check for atom_check in self.full_atom_list if (
                                    atom_check.point_type == connection.origin and atom_check.owning_cell == mirrored_target_cell)),
                                        None):
                                self.atom_count += 1
                                self.full_atom_list.append(Atom(self.atom_count,
                                                            np.add(mirrored_cell_origin, pseudo_origin_coords),
                                                            connection.origin, mirrored_target_cell))
                            mirrored_bond_list.append(
                                dict(connection=Connection(connection.partner,
                                                           connection.origin,
                                                           tuple(np.multiply(connection.cell_crossings, -1))),
                                     origin_cell=current_cell))

        for current_cell in cell_list:
            dangling_bonds = (item['connection'] for item in mirrored_bond_list if
                                  item['origin_cell'] == current_cell)
            connections = self.base_connection_list + list(dangling_bonds)
            for connection in connections:
                origin = next((atom for atom in self.full_atom_list
                               if (atom.owning_cell == current_cell
                                   and atom.point_type == connection.origin)),None)
                partner = next((atom for atom in self.full_atom_list
                                if (atom.owning_cell == tuple(np.add(current_cell, connection.cell_crossings))
                                    and atom.point_type == connection.partner)),None)
                if origin and partner:
                    origin.bonds += 1
                    partner.bonds += 1
                    self.vector_list.append(dict(origin = origin, partner = partner))
        if self.prune:
            self._prune_vectors()
        for pair in self.vector_list:
            self._populate_vector(pair)

    def _populate_vector(self, atoms):
        """
        Populate the line segment connecting two atoms with a polymer chain
        :param atoms: Dict of two Atom objects named 'origin' and 'partner';
        passed as an element of the vector list in standard usage but can be supplied directly
        """
        vector = np.subtract(atoms['partner'].position,atoms['origin'].position)
        vector_length = int(np.floor(np.sqrt((np.sum(np.square(vector)))) + self.float_safety))
        if vector_length < 1:
            vector_length = 1
        step = vector / vector_length
        for i in range(vector_length):
            if i < vector_length - 1:
                self.atom_count += 1
                new_atom = Atom(self.atom_count, np.add(atoms['origin'].position, np.multiply(step, i+1)))
                self.full_atom_list.append(new_atom)
                bond_origin = new_atom.id
            else:
                bond_origin = atoms['partner'].id
            if i == 0:
                bond_partner = atoms['origin'].id
            elif i < vector_length - 1:
                bond_partner = new_atom.id -1
            else:
                bond_partner = new_atom.id
            self.bond_count += 1
            self.full_bond_list.append(dict(bond_id=self.bond_count, point1=bond_origin, point2=bond_partner))

    def _prune_vectors(self):
        """
        Prunes the vector list by recursively removing points that only have a single connection
        """
        while True:
            dangler = next((atom for atom in self.full_atom_list if atom.bonds == 1),None)
            if not dangler:
                break
            for connection in self.vector_list:
                if connection['origin'] == dangler:
                    connection['partner'].bonds -= 1
                    self.vector_list.remove(connection)
                if connection['partner'] == dangler:
                    connection['origin'].bonds -= 1
                    self.vector_list.remove(connection)
            self.full_atom_list.remove(dangler)



