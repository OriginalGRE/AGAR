"""
Exports the Cell, Point and Connection classes for AGAR
"""
import pickle
import os
import numpy as np
from typing import Optional


class Point:
    """
    Template class for points, used in constructing Cell objects
    Cell.add_particle is recommended over direct instantiation when editing a single geometry
    Direct instantiation bypasses some useful checks, but may be helpful in scripts that create batches of geometries
    """
    __slots__ = ["point_id", "coords"]

    def __init__(self, point_id: int, x: float, y: float, z: float):
        self.point_id = point_id
        self.coords = (x, y, z)


class Connection:
    """
    Template class for connections, used in constructing Cell objects
    Cell.add_connection is recommended over direct instantiation when editing a single geometry
    Direct instantiation bypasses some useful checks, but may be helpful in scripts that create batches of geometries
    """
    __slots__ = ["origin", "partner", "cell_crossings"]

    def __init__(self, origin: int, partner: int, cell_crossings: Optional[tuple[int, int, int]] = None):
        self.origin = origin
        self.partner = partner
        self.cell_crossings = cell_crossings
        if cell_crossings is None:
            self.cell_crossings = (0, 0, 0)


class Cell:

    """
    Stores cell geometry and provides convenience methods for importing, constructing and modifying cell geometries
    """

    def __init__(self, point_list: Optional[list[Point]] = None, connection_list: Optional[list[Connection]] = None,
                 dimensions: Optional[tuple[float, float, float]] = (1, 1, 1), template: Optional[str] = None):
        """
        Constructs the cell object; either from the provided lists or from a template file.
        Other arguments are ignored if a template is provided.
        :param point_list: list of point dicts, comprised of an ID ('point') and a tuple of 3 coordinates ('position')
        :param connection_list: list of connection dicts, comprised of 2 point IDs ('origin' and 'partner') and a tuple
        representing whether the connection crosses cell edges ('cell_offset')
        :param dimensions: tuple of floats representing the dimensions of the cell in x, y, z
        :param template: filename of geometry template
        """
        if template is None:
            if point_list is None:
                self.points = []
            else:
                self.points = point_list

            if connection_list is None:
                self.connections = []
            else:
                self.connections = connection_list

            self.dims = dimensions

        else:
            self._import_cell(template)

        for point in self.points:
            if not all(np.less_equal(point.coords, self.dims)):
                print(f"WARNING: Point {point.point_id} falls outside current cell dimensions.")

    def _import_cell(self, template: str):  # Need to find a good way to import arbitrary variables
        """
        Imports cell geometry from a template file; currently set up for pickle, may support other formats in future
        :param template: filename of geometry template
        :return: None
        """
        if os.path.exists(template):
            with open(template, 'rb') as geometry:
                print("Importing geometry from {}".format(template))
                self.points, self.connections, self.dims = pickle.load(geometry)
        else:
            print("Geometry template does not exist")

    def write_cell(self, filename: str):
        """
        Dumps a pickle containing all current geometry data
        :param filename: filename for geometry file; will be overwritten if filename already exists
        :return: None
        """
        with open(filename, 'wb') as geometry:
            pickle.dump((self.points, self.connections, self.dims), geometry)

    def number_of_points(self) -> int:
        """
        Counts the number of points in the cell instance
        :return: int: number of points
        """
        return len(self.points)

    def number_of_connections(self) -> int:
        """
        Counts the number of connections in the cell instance
        :return: int: number of connections
        """
        return len(self.connections)

    def set_cell_dims(self, x_new: float, y_new: float, z_new: float, rescale: Optional[bool] = False):
        """
        Change the dimensions of the cell
        :param x_new: new cell dimension in x
        :param y_new: new cell dimension in y
        :param z_new: new cell dimension in z
        :param rescale: bool indicating whether to scale existing particle coordinates along with cell dimensions
        :return: None
        """
        new_dims = (x_new, y_new, z_new)
        if rescale:
            scaling_factors = np.divide(new_dims, self.dims)
            print("Rescaling particle positions")

            for point in self.points:
                point.coords = np.multiply(point.coords, scaling_factors)

        self.dims = new_dims

    def add_point(self, x: float, y: float, z: float, point_id: Optional[int] = None):
        """
        Adds a point to the cell by instantiating a Point object and appending it to the point list.
        Prints warnings if the new point falls outside the cell dimensions
        :param x: coordinate along dimension x
        :param y: coordinate along dimension y
        :param z: coordinate along dimension z
        :param point_id: Point identifier for bonding purposes; defaults to the number of points
        (if point_id is never explicitly set, point ids will be contiguous integers in order of creation)
        :return:
        """
        if point_id is None:
            point_id = self.number_of_points() + 1
        self.points.append(Point(point_id, x, y, z))
        if not all(np.less_equal((x, y, z), self.dims)):
            print(f"Warning: New point {point_id} falls outside current cell dimensions.")

    def add_connection(self, origin_id: int, partner_id: int, cell_offset: Optional[tuple[int, int, int]] = (0, 0, 0)):
        """
        Adds a new connection to the cell by instantiating a Connection object and adding it to the connection list.
        A warning is printed if the connection includes nonexistent particles
        :param origin_id: ID of the first point in the connection
        :param partner_id: ID of the second point in the connection
        :param cell_offset: Tuple of ints; defines whether the connection crosses cell boundaries
        Defaults to (0, 0, 0), i.e. a connection within the cell
        :return: None
        """
        self.connections.append(Connection(origin_id, partner_id, cell_offset))
        origin_exists = (next((point for point in self.points if point.point_id == origin_id), None) is not None)
        partner_exists = (next((point for point in self.points if point.point_id == partner_id), None) is not None)

        if not (origin_exists and partner_exists):
            print("Warning: at least one connected point does not currently exist")
