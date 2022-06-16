import sys
import os
sys.path.append(os.path.dirname(os.path.realpath(__file__)))
from readlayer import layer_trj, layer_vtk, layer_data, singleAggregate


class particlelayer:

    def __init__(self, inputfile, nodescription=True):
        self.file = inputfile
        self._format = None
        self._layer = None
        self._lib = {'trj': layer_trj, 'vtk': layer_vtk, 'datafile': layer_data, 'singleAggregate': singleAggregate}
        if not nodescription:
            self.printInfo()

    @property
    def layer(self):
        if self._layer is None:
            try:
                self._layer = self._lib[self.format](self.file)
            except KeyError:
                print("Could not find key {} for file {}".format(self.format, self.file))
                sys.exit()
        return self._layer

    @property
    def format(self):
        if self._format is None:
            self.getFormat()
        return self._format

    def getFormat(self):
        if self.file.endswith('vtk'):
            with open(self.file) as f:
                lines = f.readlines()
            if 'vtk' not in lines[0]:
                raise IOError('Cannot determine fileformat of %s (maybe vtk?)' % self.file)
            self._format = 'vtk'

        elif self.file.endswith('trj') or 'dump' in self.file:
            with open(self.file) as f:
                lines = f.readlines()
            if 'ITEM:' not in lines[0]:
                raise IOError('Cannot determine fileformat of %s (maybe dump file?)' % self.file)
            self._format = 'trj'
        elif 'data' in self.file:
            with open(self.file) as f:
                lines = f.readlines()
            if 'data' not in lines[0]:
                raise IOError('Cannot determine fileformat of %s (maybe datafile?)' % self.file)
            self._format = 'datafile'

        elif self.file.endswith('.dat'):
            with open(self.file) as f:
                lines = f.readlines()
            if lines[0].startswith('# Aggregate data file'):
                self._format = 'singleAggregate'

    def printInfo(self):
        print('\n### Reading Particle Structure File')
        print('# Filename: {}'.format(self.file))
        string = ''
        string += ('# File Format: {}'.format(self.format))
        string += ('\n# Units: {}'.format(self.layer.units))
        string += ('\n# Primary Particles: {:,}'.format(self.layer.particles))
        string += ('\n# Aggregates: {:d}'.format(int(max(self.layer.aggregates))))
        string += ('\n# Box Width: {:1.3e} x {:1.3e}'.format(self.layer.box[0, 1] - self.layer.box[0, 0], self.layer.box[1, 1] - self.layer.box[1, 0]))
        string += ('\n# Film Height: {:1.3e}'.format(self.layer.box[2, 1] - self.layer.box[2, 0]))
        string += '\n### Done\n\n'

        print(string)
