from Cell import Cell
from Lattice_experimental import Lattice
from LammpsDataWriter import LammpsDataWriter
import numpy as np

small_lattice = Lattice(Cell(template='diamond.pickle'), (2,2,2), dangling=True, symmetric=True, prune=True, scale=32/np.sqrt(3))
LammpsDataWriter(small_lattice, 'diamond_gel_222.data')

large_lattice = Lattice(Cell(template='diamond.pickle'), (3,3,3), dangling=True, symmetric=True, prune=True, scale=32/np.sqrt(3))
LammpsDataWriter(large_lattice, 'diamond_gel_333.data')
