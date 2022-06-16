import sys
import os
from layer_main import layer
import numpy as np


class layer_trj(layer):

    def __init__(self, inputfile, filetype='dda'):
        layer.__init__(self)
        self.file = inputfile
        self._lines = None
        self._lenheader = None
        self._units = None
        self._filetype = filetype
        self._lib = None
        self._particles = None
        self._datalib = None
        self._rmax = None
        self._data = None

        self.getFileProperties()

    @property
    def lines(self):
        if self._lines is None:
            with open(self.file, 'r') as f:
                self._lines = f.readlines()
        return self._lines

    def getBox(self):
        """
        Reads the box dimensions
        :return: The Box dimensions in [-x, +x, -y, +y, -z, +z]
        """
        for line in range(len(self.lines)):
            if self.lines[line].startswith("ITEM: BOX"):
                break
        self._box = np.zeros((3, 2))
        self._box[0, 0] = float(self.lines[line + 1].split()[0])
        self._box[0, 1] = float(self.lines[line + 1].split()[1])
        self._box[1, 0] = float(self.lines[line + 2].split()[0])
        self._box[1, 1] = float(self.lines[line + 2].split()[1])
        self._box[2, 0] = float(self.lines[line + 3].split()[0])
        self._box[2, 1] = float(self.lines[line + 3].split()[1])
        if self._units == 'si':
            self._box *= 1e9

    @property
    def headlines(self):
        """
        Get the length of the headerfile to only extract the data
        :return: The amount of header lines in the file
        """
        if self._lenheader is None:
            for line in range(len(self.lines)):
                if self.lines[line].startswith('ITEM: ATOMS'):
                    self._lenheader = line + 1
                    break
        return self._lenheader

    def getFileProperties(self, verbose=True):
        if self.box[0, 1] - self.box[0, 0] > 1:
            self._units = 'nano'
        else:
            self._units = 'si'
        if 'omega' in self.lines[self.headlines - 1]:
            self._filetype = 'liggghts'
        else:
            self._filetype = 'dda'
        print('    File {} using units {} and type {}'.format(self.file, self._units, self._filetype))
        if self._units == 'si':
            print('    Converting to nano')
            self.data[:, self.lib['x']] *= 1e9
            self.data[:, self.lib['y']] *= 1e9
            self.data[:, self.lib['z']] *= 1e9
            self.data[:, self.lib['r']] *= 1e9
            self._box *= 1e9
            self._units = 'nano'

    @property
    def filetype(self):
        return self._filetype

    @filetype.setter
    def filetype(self, filetype):
        if filetype not in ('liggghts', 'dda'):
            raise ValueError('Only filetypes liggghts and dda are supported')
        self._filetype = filetype

    @property
    def lib(self):
        if self._lib is None:
            if self.filetype == 'dda':
                self._lib = {"id": 0, "aggregate": 1, "cluster": 2, "type": 3, "x": 4, "y": 5, "z": 6,
                             "r": 7, "c_sa": 8, "c_da": 9, "c_tot": 10, "c_neigh": 11, "c_same": 12, "brokenBondFlag": 13, "percolationPath": 14, "relDisplacementLast": 15, "relDisplacementTotal": 16}

            elif self.filetype == 'liggghts':
                self._lib = {"id": 0, "aggregate": 1, "type": 2, "r": 3, "x": 4, "y": 5, "z": 6, "ix": 7,
                             "iy": 8, "iz": 9, "vx": 10, "vy": 11, "vz": 12, "fx": 13, "fy": 14, "fz": 15,
                             "fmag": 16, "omegax": 17, "omegay": 18, "omegaz": 19, "c_sa": 20, "c_da": 21,
                             "c_tot": 22, "brokenBondFlag": 23, "percolationPath": 24, "relDisplacementLast": 25, "relDisplacementTotal": 26}
        return self._lib

    @property
    def particles(self):
        """
            :return: The amount of particles within the box
        """
        if self._particles is None:
            self._particles = len(self.lines) - self.headlines
        return self._particles

    @particles.setter
    def particles(self, amount):
        if type(amount) == int:
            self._particles = amount
        else:
            raise ValueError('Only integers can be used for amount of particles')

    @property
    def datalib(self):
        if self._datalib is None:
            lib = self.lines[self.headlines - 1].split()[2:]
            for dat in range(len(lib)):
                if lib[dat] == 'mol':
                    lib[dat] = 'aggregate'
                if self.filetype == 'dda' and lib[dat] == 'type' and  \
                        'aggregate' not in lib:
                    lib[dat] = 'aggregate'
                if lib[dat].lower() == 'radius':
                    lib[dat] = 'r'
            self._datalib = {}
            for dat in range(len(lib)):
                self._datalib[lib[dat]] = dat
        return self._datalib

    def getData(self):
        if self._data is None:
            self._data = np.zeros((self.particles, len(self.lib)))
            for lin in range(len(self.lines) - self.headlines):
                tmp = self.lines[lin + self.headlines].split()
                for y in range(len(tmp)):
                    self._data[lin, self.lib[self.getKey(y)]] = float(tmp[y])
            del self._lines
            self._lines = None

    def getKey(self, entry):
        for key in self.datalib:
            if entry == self.datalib[key]:
                return key


class layer_vtk(layer):

    def __init__(self, inputfile):
        layer.__init__(self)
        self.file = inputfile
        self._lines = None
        self._lenheader = None
        self._units = None
        self._filetype = None
        self._lib = None
        self._particles = None
        self._datalib = None
        self._rmax = None
        self._data = None

        self.getFileProperties()

    @property
    def lines(self):
        if self._lines is None:
            with open(self.file, 'r') as f:
                self._lines = f.readlines()
        return self._lines

    def getBox(self):
        """
        Reads the box dimensions
        :return: The Box dimensions in [-x, +x, -y, +y, -z, +z]
        """
        boxfile = self.file.split('.')[0] + '_box.vtk'
        if not os.path.exists(boxfile):
            boxfile = self.file.split('.')[0] + '_Box.vtk'
        try:
            with open(boxfile) as f:
                lines = f.readlines()
        except FileNotFoundError:
            raise FileNotFoundError('File could not be found. Please use inputfile_box.vtk')
        self._box = np.zeros((3, 2))
        for line in range(len(lines)):
            if lines[line].startswith("X_COORDINATES"):
                self._box[0, 0] = float(lines[line + 1].split()[0])
                self._box[0, 1] = float(lines[line + 1].split()[1])
            elif lines[line].startswith("Y_COORDINATES"):
                self._box[1, 0] = float(lines[line + 1].split()[0])
                self._box[1, 1] = float(lines[line + 1].split()[1])
            elif lines[line].startswith("Z_COORDINATES"):
                self._box[2, 0] = float(lines[line + 1].split()[0])
                self._box[2, 1] = float(lines[line + 1].split()[1])
                break

    def getFileProperties(self):
        for lin in range(len(self.lines)):
            if self.lines[lin].startswith('POINTS'):
                self._particles = int(self.lines[lin].split()[1])
                if float(self.lines[lin + 1].split()[0]) < 1e-3:
                    self._units = 'si'
                else:
                    self._units = 'nano'
                break
        self._lib = {'x': 0, 'y': 1, 'z': 2}
        count = 3
        for lin in self.lines:
            if lin.startswith('SCALARS'):
                self._lib[lin.split()[1]] = count
                count += 1
        if 'brokenBondFlag' not in self._lib.keys():
            self._lib['brokenBondFlag'] = count
            count += 1
        if 'percolationPath' not in self._lib.keys():
            self._lib['percolationPath'] = count
            count += 1
        if 'c_sa' not in self._lib.keys():
            self._lib['c_sa'] = count
            count += 1
            self._lib['c_da'] = count
            count += 1
            self._lib['c_tot'] = count
            count += 1
        print('File using units {}'.format(self._units))
        if self._units == 'si':
            print('Converting to nano')
            self.data[:, self.lib['x']] *= 1e9
            self.data[:, self.lib['y']] *= 1e9
            self.data[:, self.lib['z']] *= 1e9
            self.data[:, self.lib['r']] *= 1e9

    @property
    def lib(self):
        return self._lib

    def getData(self):
        self._data = np.zeros((self.particles, len(self.lib)))
        for lin in range(len(self.lines)):
            if self.lines[lin].startswith('POINTS'):
                start = lin + 1
                for pp in range(self.particles):
                    self._data[pp, self.lib['x']] = float(self.lines[pp + start].split()[0])
                    self._data[pp, self.lib['y']] = float(self.lines[pp + start].split()[1])
                    self._data[pp, self.lib['z']] = float(self.lines[pp + start].split()[2])
            elif self.lines[lin].startswith('SCALARS'):
                start = lin + 2
                prop = self.lines[lin].split()[1]
                for pp in range(self.particles):
                    try:
                        self._data[pp, self.lib[prop]] = float(self.lines[pp + start].split()[0])
                    except IndexError:
                        print('Error at property {}'.format(prop))
                        print('Should be {} values, cannot access value {}'.format(self.particles, pp + 1))
                        print('List of properties')
                        print(self.lib)
                        raise IndexError
        # Check if any index value has gotten wrong
        minID = min(self._data[:, self.lib['id']])
        mintype = min(self._data[:, self.lib['type']])
        minaggregate = min(self._data[:, self.lib['aggregate']])
        for pp in range(len(self._data)):
            self._data[pp, self.lib['id']] = self._data[pp, self.lib['id']] - minID + 1
            self._data[pp, self.lib['type']] = self._data[pp, self.lib['type']] - mintype + 1
            self._data[pp, self.lib['aggregate']] = self._data[pp, self.lib['aggregate']] - minaggregate + 1


class layer_data(layer):

    def __init__(self, inputfile):

        layer.__init__(self)
        self.file = inputfile
        self._lines = None
        self._units = None
        self._format = None
        self._lib = None
        self._datalib = None

        self.getFileProperties()

    @property
    def lines(self):
        if self._lines is None:
            with open(self.file) as f:
                self._lines = f.readlines()
        return self._lines

    def getFileProperties(self):
        for lin in self.lines:
            if "atoms" in lin:
                self._particles = int(lin.split()[0])
                break
        if "si" in self.lines[0]:
            self._units = 'si'
        elif "nano" in self.lines[0]:
            self._units = 'nano'
        if "hybrid granular bond/gran" in self.lines[0]:
            self._format = "hybrid granular bond/gran"
        elif "hybrid granular molecular" in self.lines[0]:
            self._format = "hybrid granular molecular"
        print('File using units {}'.format(self._units))
        if self._units == 'si':
            print('Converting to nano')
            self.data[:, self.lib['x']] *= 1e9
            self.data[:, self.lib['y']] *= 1e9
            self.data[:, self.lib['z']] *= 1e9
            self.data[:, self.lib['r']] *= 1e9

    def getBox(self):
        for lin in range(len(self.lines)):
            if "xlo" in self.lines[lin]:
                break
        self._box = np.zeros((3, 2))
        self._box[0, 0] = float(self.lines[lin + 0].split()[0])
        self._box[0, 1] = float(self.lines[lin + 0].split()[1])
        self._box[1, 0] = float(self.lines[lin + 1].split()[0])
        self._box[1, 1] = float(self.lines[lin + 1].split()[1])
        self._box[2, 0] = float(self.lines[lin + 2].split()[0])
        self._box[2, 1] = float(self.lines[lin + 2].split()[1])

    @property
    def datalib(self):
        if self._datalib is None:
            self._datalib = {'id': 0, 'type': 1, 'x': 2, 'y': 3, 'z': 4, 'r': 5, 'aggregate': 7}
        return self._datalib

    @property
    def lib(self):
        if self._lib is None:
            self._lib = {'id': 0, 'type': 1, 'x': 2, 'y': 3, 'z': 4, 'r': 5, 'aggregate': 6, 'cluster': 7, 'brokenBondFlag': 8, "relDisplacementLast": 9, "relDisplacementTotal": 10}
        return self._lib

    def getData(self):
        start = -1
        for lin in range(len(self.lines)):
            if self.lines[lin].startswith('Atoms'):
                if self.lines[lin + 1] == '\n':
                    start = lin + 2
                else:
                    start = lin + 1
        self._data = np.zeros((self.particles, len(self.lib)))
        for lin in range(len(self.lines) - start):
            tmp = self.lines[lin + start].split()
            for entry in self.datalib:
                self._data[lin, self.lib[entry]] = float(tmp[self.datalib[entry]])
        # Diameter is given in datafile
        for pp in range(len(self._data)):
            self._data[pp, self.lib['r']] = (self._data[pp, self.lib['r']] - 0.15e-9) / 2


class singleAggregate(layer):

    def __init__(self, inputfile):
        layer.__init__(self)
        self.file = inputfile
        self._lines = None
        self._lenheader = None
        self._units = None
        self._lib = None
        self._particles = None
        self._data = None
        self._ppsizedistribution = None
        self._ppsizedistributionfile = None
        self._D_f = None
        self._k_f = None
        self._box = None
        self.getFileProperties()

    @property
    def box(self):
        if self._box is None:
            self._box = np.zeros((3, 2))
            self._box[0, 0] = min(self.data[:, self.lib['x']])
            self._box[0, 1] = max(self.data[:, self.lib['x']])
            self._box[1, 0] = min(self.data[:, self.lib['y']])
            self._box[1, 1] = max(self.data[:, self.lib['y']])
            self._box[2, 0] = min(self.data[:, self.lib['z']])
            self._box[2, 1] = max(self.data[:, self.lib['z']])
        return self._box

    @property
    def lines(self):
        if self._lines is None:
            with open(self.file, 'r') as f:
                self._lines = f.readlines()
        return self._lines

    @property
    def headlines(self):
        """
        Get the length of the headerfile to only extract the data
        :return: The amount of header lines in the file
        """
        if self._lenheader is None:
            for line in range(len(self.lines)):
                if self.lines[line].startswith('PARTICLE DATA'):
                    self._lenheader = line + 1
                    break
        return self._lenheader

    def getData(self):
        if self._data is None:
            self._data = np.zeros((self.particles, len(self.lib)))
            for lin in range(len(self.lines) - self.headlines):
                tmp = self.lines[lin + self.headlines].split()
                for y in range(len(tmp)):
                    self._data[lin, self.lib[self.getKey(y)]] = float(tmp[y])
                if self.lines[lin + self.headlines + 1].startswith('END'):
                    break
            if sum(self.data[:, self.lib['aggregate']]) == 0:
                self.data[:, self.lib['aggregate']] = 1
            if sum(self.data[:, self.lib['type']]) == 0:
                self.data[:, self.lib['type']] = 1

    def getKey(self, entry):
        for key in self.lib:
            if entry == self.lib[key]:
                return key

    @property
    def particles(self):
        """
            :return: The amount of particles within the box
        """
        if self._particles is None:
            self._particles = int(self.lines[1].split(':')[-1])
        return self._particles

    @particles.setter
    def particles(self, amount):
        if type(amount) == int:
            self._particles = amount
        else:
            raise ValueError('Only integers can be used for amount of particles')

    @property
    def lib(self):
        if self._lib is None:
            self._lib = {"id": 0, "cluster": 1, "x": 2, "y": 3, "z": 4,
                         "r": 5, "aggregate": 6, "type": 7}
        return self._lib

    def getFileProperties(self):
        if self.data[0, self.lib['r']] < 1e-3:
            self._units = 'si'
        else:
            self._units = 'nano'
        for lin in range(self.headlines):
            tmp = self.lines[lin].split(':')
            if tmp[0] == 'particles':
                self._particles = int(tmp[1])
            elif tmp[0] == 'pp size distribution':
                self._ppsizedistribution = tmp[1]
            elif tmp[0] == 'pp size distribution file':
                self._ppsizedistributionfile = tmp[1]
            elif tmp[0] == 'fractal dimension':
                self._D_f = float(tmp[1])
            elif tmp[0] == 'fractal prefactor':
                self._k_f = float(tmp[1])

    @property
    def ppsizedistribution(self):
        return self._ppsizedistribution

    @property
    def ppsizedistributionfile(self):
        return self._ppsizedistributionfile

    @property
    def D_f(self):
        return self._D_f

    @property
    def k_f(self):
        return self._k_f
