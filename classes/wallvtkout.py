import sys
import os
import numpy as np
sys.path.append(os.path.dirname(os.path.realpath(__file__)))
from convert_units import conversion


class wallvtkout:
    def __init__(self, inputfile, units='si'):

        self.file = inputfile
        with open(self.file) as f:
            self.lines = f.readlines()

        self._points = None
        self._pointcoordinates = None
        self._timestep = None
        self._walls = None
        self._cells = None
        self._cellindices = None
        self._celltypes = None
        self._celltype = None
        self._pressure = None
        self._shearstress = None
        self._planeXCoordinate = None
        self._planeYCoordinate = None
        self._planeZCoordinate = None
        self._boxVolume = None
        self.convert = conversion(unit_from=units, unit_to='si').convert

    @property
    def timestep(self):
        """
        Return the timestep of the current file
        """
        if self._timestep is None:
            self._timestep = int(self.file.split('_')[-1].split('.')[0])
        return self._timestep

    @property
    def points(self):
        """
        Return the amount of vertices in the grid
        """
        if self._points is None:
            for lin in self.lines:
                if lin.startswith('POINTS'):
                    self._points = int(lin.split()[1])
                    break
        return self._points

    @property
    def pointcoordinates(self):
        """
        Return the coordinates of the vertices
        """
        if self._pointcoordinates is None:
            for lin in range(len(self.lines)):
                if self.lines[lin].startswith('POINTS'):
                    start = lin + 1
                    break
            self._pointcoordinates = np.zeros((self.points, 3))
            for coord in range(self.points):
                self._pointcoordinates[coord] = self.lines[coord + start].split()
        return self._pointcoordinates

    @property
    def walls(self):
        """
        Return the amount of different walls included in this file
        """
        if self._walls is None:
            self._walls = int(len(set(self.pointcoordinates[:, 2])) / 2)
        return self._walls

    @property
    def cells(self):
        """
        Return the amount of cells included in this file
        """
        if self._cells is None:
            for lin in self.lines:
                if lin.startswith('CELLS'):
                    self._cells = int(lin.split()[1])
                    break
        return self._cells

    @property
    def cellindices(self):
        """
        Return the indices for the cells
        """
        if self._cellindices is None:
            for lin in range(len(self.lines)):
                if self.lines[lin].startswith('CELLS'):
                    start = lin + 1
                    break
            self._cellindices = np.zeros((self.cells, 4))
            for cd in range(self.cells):
                self._cellindices[cd] = self.lines[cd + start].split()
        return self._cellindices

    @property
    def celltypes(self):
        """
        Return the different cell types
        """
        if self._celltypes is None:
            for lin in range(len(self.lines)):
                if self.lines[lin].startswith('CELL_TYPES'):
                    start = lin + 1
                    break
            self._celltypes = np.zeros(self.cells)
            for ty in range(self.cells):
                self._celltypes[ty] = int(self.lines[ty + start].split()[0])
        return self._celltypes

    @property
    def celldata(self):
        return {'pressure': self.pressure, 'shearstress': self.shearstress}

    @property
    def pressure(self):
        """
        Return the pressure values per cell in the format:
        pressure[x,y] = pressure of cell y in wall x
        """
        if self._pressure is None:
            for lin in range(len(self.lines)):
                if self.lines[lin].startswith('SCALARS normal_stress') or \
                   self.lines[lin].startswith('SCALARS pressure'):
                    start = lin + 2
                    break
            self._pressure = np.zeros((self.walls, int(self.cells / self.walls)))
            for data in range(int(self.cells / self.walls)):
                self._pressure[0, data] = float(self.lines[data + start].split()[0]) * self.convert('pressure')
                if self.walls > 1:
                    self._pressure[1, data] = float(self.lines[data + start + int(self.cells / self.walls)].split()[0]) * self.convert('pressure')
        return self._pressure

    @property
    def shearstress(self):
        """
        Return the shearsress values per cell in the format:
        shearstress[x,y] = shearstress of cell y in wall x
        """
        if self._shearstress is None:
            for lin in range(len(self.lines)):
                if self.lines[lin].startswith('SCALARS shear_stress') or \
                   self.lines[lin].startswith('SCALARS shearstress'):
                    start = lin + 2
                    break
            self._shearstress = np.zeros((self.walls, int(self.cells / self.walls)))
            for data in range(int(self.cells / self.walls)):
                self._shearstress[0, data] = self.lines[data + start].split()[0]
                if self.walls > 1:
                    self._shearstress[1, data] = self.lines[data + start + int(self.cells / self.walls)].split()[0]
        return self._shearstress

    @property
    def meanvalue(self):
        """
        Return the mean value of a property for all cells with magnitude > 0
        """
        meanpressure = np.zeros(self.walls)
        meanstress = np.zeros(self.walls)

        for wa in range(self.walls):
            tmppressure = []
            tmpstress = []
            for dat in range(int(self.cells / self.walls)):
                if self.celldata['pressure'][wa, dat] > 0:
                    tmppressure.append(self.celldata['pressure'][wa, dat])
                    tmpstress.append(self.celldata['shearstress'][wa, dat])
            if len(tmppressure) > 0:
                meanpressure[wa] = np.mean(tmppressure)
                meanstress[wa] = np.mean(tmpstress)
            else:
                meanpressure[wa] = 0
                meanstress[wa] = 0

        return {'pressure': meanpressure, 'shearstress': meanstress}

    """
    The coordinates of the inner planes of the walls (planeZCoordinate) as well as the outer boundaries
    are necessary to know to account for the volume between the walls
    """

    @property
    def planeZCoordinate(self):
        if self._planeZCoordinate is None:
            self._planeZCoordinate = np.array(sorted(list(set(self.pointcoordinates[:, 2]))))
        return self._planeZCoordinate

    @property
    def planeXCoordinate(self):
        if self._planeXCoordinate is None:
            self._planeXCoordinate = np.array([min(self.pointcoordinates[:, 0]), max(self.pointcoordinates[:, 0])])
        return self._planeXCoordinate

    @property
    def planeYCoordinate(self):
        if self._planeYCoordinate is None:
            self._planeYCoordinate = np.array([min(self.pointcoordinates[:, 1]), max(self.pointcoordinates[:, 1])])
        return self._planeYCoordinate

    @property
    def boxVolume(self):
        if self._boxVolume is None:
            self._boxVolume = (self.planeXCoordinate[1] - self.planeXCoordinate[0]) *  \
                              (self.planeYCoordinate[1] - self.planeYCoordinate[0]) * \
                              (self.planeZCoordinate[2] - self.planeZCoordinate[1])
        return self._boxVolume
