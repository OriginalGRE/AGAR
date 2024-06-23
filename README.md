# AGAR

AGAR is a set of Python classes for creating and modifying input files for molecular dynamics simulations, particularly simulations of polymer networks. The current core classes are geared towards Kremer-Grest simulations in LAMMPS, but written with the goal of being easy to adapt and extend to other systems and software. AGAR refers to agar gel, and is an acronym of Automated Gel Architecture Replicator.

## Overview

AGAR produces a system configuration for molecular dynamics simulations in three steps. The classes in `Cell.py` are used to define a volume of space, the network crosslinks contained in that volume, and their connectivity. The `Lattice` class is used to scale a `Cell` geometry to the desired size, and tile copies of the `Cell` together into a larger structure. `Lattice` also provides methods to create polymer chains along connections, and options for handling "dangling" connections. Finally, `LammpsDataWriter` reads the atom and bond lists defined by a `Lattice` and writes out a corresponding LAMMPS data file. 

### Cell

`Cell.py` contains the `Point`, `Connection` and `Cell` classes, with `Point` and `Connection` primarily serving as a standard structure for geometry input. The properties of a `Point` are a numerical `point_id` which is used in constructing connectivity; and `coords` which define its position in the cell. The properties of a `Connection` are `origin` and `partner`, which correspond to the `point_id` of the points it connects; and `cell_crossings`, which defines how many times the connection crosses a cell boundary in each dimension. The `Cell` class defines the geometry of a volume of polymer gel. The properties of a `Cell` object are `dims`, which defines the dimensions of the (orthogonal) cell in x, y and z; `points`, which is a list of `Point` objects that define crosslinks contained in the cell; and `connections`, which is (predictably) a list of `Connection` objects that describe the connectivity between points in the network. 

`Cell` provides various convenience methods for defining structures. `write_cell()` creates a pickle file containing the properties of the `Cell`, which can be imported by using the `template` keyword in the constructor. `add_point()` creates a `Point` object with the specified coordinates and appends it to the point list. `add_connection()` likewise creates a `Connection` between two points and appends it to the connection list. `set_cell_dims()` alters the dimensions of the cell, with the option to rescale `Point` coordinates along with it. 

`Point` and `Connection` objects can be instantiated independent of a `Cell`; this may be convenient for automating the creation of batches of different geometries. Similarly, point and connection lists can be passed to the `Cell` constructor as arguments.

**At least one connection outside the cell is required to produce a contiguous network!**

### Lattice

`Lattice.py` contains the `Atom` and `Lattice` classes, with `Atom` serving as a data structure for internal use by `Lattice`; each `Atom` describes a particle in the final data file. The `Lattice` class generates a full gel structure from a `Cell` object according to the `size` argument, which is a tuple of 3 integers describing how many times to tile the `Cell` in each dimension. For instance, a `size` of (1, 1, 3) would correspond to three cells in a row, stacked in z. The optional `scale` argument can be used to scale the dimensions and point positions of the `Cell`. The optional `dangling`, `symmetric`, and `prune` arguments define the handling of various edge cases. `dangling` defines whether connections that originate within the `Lattice` but connect to a (theoretical) particle outside it are created or not; it defaults to `False`. If `dangling` is allowed, `symmetric` defines whether dangling connections that only occur on one side of the `Lattice` are mirrored on the other; it defaults to `True`, but this has no effect unless `dangling` is enabled. If `prune` is `True`, after the full network structure is generated, points with only one connection are recursively removed by the `_prune_vectors()` method, leaving a structure with no free chain ends; `prune` defaults to `False`. 

The main functionality of the `Lattice` class is in the `build_lists()` method. This method loops over all cells to construct a list of relevant points and connections, then creates polymer chains along all connections by calling the _populate_vector() method, and constructs a list of atoms and bonds to be written out in the final data file.

### LAMMPSDataWriter

`LammpsDataWriter` is a fairly simple class that accepts a `Lattice` object, reads its `full_atom_list` and `full_bond_list` properties, and outputs a corresponding LAMMPS data file. 
