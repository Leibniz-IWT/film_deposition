class writetrj:

    def __init__(self, obj, units='si'):

        self.obj = obj
        self.units = units

        if self.units == self.obj.units:
            self.scaleup = 1
        elif self.units == 'si' and self.obj.units == 'nano':
            self.scaleup = 1e-9
        elif self.units == 'nano' and self.obj.units == 'si':
            self.scaleup = 1e9
        else:
            self.scaleup = 1
        self._propOut = None
        self._outlib = None
        self._scale = None

    @property
    def outlib(self):
        if self._outlib is None:
            self._outlib = {'id': '%i',
                            'type': '%i',
                            'aggregate': '%i',
                            'cluster': '%i',
                            'x': '%e',
                            'y': '%e',
                            'z': '%e',
                            'r': '%e',
                            'brokenBondFlag': '%i',
                            'percolationPath': '%i',
                            'c_sa': '%i',
                            'c_da': '%i',
                            'c_tot': '%i',
                            "relDisplacementLast": '%1.3f',
                            "relDisplacementTotal": '%1.3f',
                            }
        return self._outlib

    @outlib.setter
    def outlib(self, key, value):
        self._outlib[key] = value

    @property
    def propOut(self):
        if self._propOut is None:
            self._propOut = {'id': 'id',
                             'type': 'type',
                             'aggregate': 'aggregate',
                             'cluster': 'cluster',
                             'x': 'x',
                             'y': 'y',
                             'z': 'z',
                             'r': 'Radius',
                             'brokenBondFlag': 'brokenBondFlag',
                             'percolationPath': 'percolationPath',
                             'c_sa': 'c_sa',
                             'c_da': 'c_da',
                             'c_tot': 'c_tot',
                             "relDisplacementLast": "relDisplacementLast",
                             "relDisplacementTotal": "relDisplacementTotal",
                             }
        return self._propOut

    @propOut.setter
    def propOut(self, key, value):
        self._propOut[key] = value

    @property
    def scale(self):
        if self._scale is None:
            self._scale = {'id': lambda x: x,
                           'type': lambda x: x,
                           'aggregate': lambda x: x,
                           'cluster': lambda x: x,
                           'x': lambda x: x * self.scaleup,
                           'y': lambda x: x * self.scaleup,
                           'z': lambda x: x * self.scaleup,
                           'r': lambda x: x * self.scaleup,
                           'brokenBondFlag': lambda x: x,
                           'percolationPath': lambda x: x,
                           'c_sa': lambda x: x,
                           'c_da': lambda x: x,
                           'c_tot': lambda x: x,
                           "relDisplacementLast": lambda x: x,
                           "relDisplacementTotal": lambda x: x,
                           }
        return self._scale

    @scale.setter
    def scale(self, key, value):
        self._scale[key] = value

    def write(self, outputfile, addlib=None, custom_props=None, **kwargs):
        if not addlib is None:
            for lib in addlib:
                self.propOut[lib[0]] = lib[0]
                self.outlib[lib[0]] = lib[1]
                self.scale[lib[0]] = lambda x: x * lib[2]
        if custom_props is None:
            props = []
            prop = ['id', 'type', 'aggregate', 'cluster', 'x', 'y', 'z', 'r', 'brokenBondFlag', 'percolationPath', 'c_sa', 'c_da', 'c_tot', "relDisplacementLast", "relDisplacementTotal"]
            for p in prop:
                if p in self.obj.lib:
                    props.append(p)
        else:
            props = []
            for p in custom_props:
                if p in self.outlib:
                    props.append(p)

        string = 'ITEM: TIMESTEP\n0\n'
        string += 'ITEM: NUMBER OF ATOMS\n%i\n' % (self.obj.particles)
        string += 'ITEM: BOX BOUNDS pp pp pp\n'
        string += '%1.3e %1.3e\n' % (self.scale['x'](self.obj.box[0, 0]), self.scale['x'](self.obj.box[0, 1]))
        string += '%1.3e %1.3e\n' % (self.scale['x'](self.obj.box[1, 0]), self.scale['x'](self.obj.box[1, 1]))
        string += '%1.3e %1.3e\n' % (self.scale['x'](self.obj.box[2, 0]), self.scale['x'](self.obj.box[2, 1]))
        string += 'ITEM: ATOMS '
        for prop in props:
            string += self.propOut[prop] + ' '
        string += '\n'
        for pp in range(self.obj.particles):
            for prop in props:
                string += self.outlib[prop] % self.scale[prop](self.obj.data[pp, self.obj.lib[prop]])
                string += ' '
            string += '\n'

        with open(outputfile, 'w') as f:
            f.write(string)


class writeVTK:

    def __init__(self, obj, units='si'):

        self.obj = obj
        self.units = units

        if self.units == self.obj.units:
            self.scaleup = 1
        elif self.units == 'si' and self.obj.units == 'nano':
            self.scaleup = 1e-9
        elif self.units == 'nano' and self.obj.units == 'si':
            self.scaleup = 1e9
        self._outlib = None
        self._scale = None

    @property
    def scale(self):
        if self._scale is None:
            self._scale = {'id': lambda x: x + 1,
                           'type': lambda x: x + 1,
                           'aggregate': lambda x: x + 1,
                           'cluster': lambda x: x + 1,
                           'x': lambda x: x * self.scaleup,
                           'y': lambda x: x * self.scaleup,
                           'z': lambda x: x * self.scaleup,
                           'r': lambda x: x * self.scaleup,
                           'brokenBondFlag': lambda x: x,
                           'percolationPath': lambda x: x,
                           'c_sa': lambda x: x,
                           'c_da': lambda x: x,
                           'c_tot': lambda x: x,
                           "relDisplacementLast": lambda x: x,
                           "relDisplacementTotal": lambda x: x,
                           }
        return self._scale

    @scale.setter
    def scale(self, key, value):
        self._scale[key] = value

    @property
    def outlib(self):
        if self._outlib is None:
            self._outlib = {'id': '%1.1f',
                            'type': '%1.1f',
                            'aggregate': '%1.1f',
                            'cluster': '%1.1f',
                            'x': '%e',
                            'y': '%e',
                            'z': '%e',
                            'r': '%e',
                            'brokenBondFlag': '%1.1f',
                            'percolationPath': '%i',
                            'c_sa': '%i',
                            'c_da': '%i',
                            'c_tot': '%i',
                            "relDisplacementLast": '%1.3f',
                            "relDisplacementTotal": '%1.3f',
                            }
        return self._outlib

    def header(self):
        string = '# vtk DataFile Version 2.0\n'
        string += 'Generated by dda\n'
        string += 'ASCII\n'
        string += 'DATASET POLYDATA\n'
        return string

    def coordinates(self):
        string = 'POINTS %i float\n' % self.obj.particles
        for particle in range(self.obj.particles):
            string += '%e %e %e\n' % (self.scale['x'](self.obj.data[particle, self.obj.lib['x']]),
                                      self.scale['x'](self.obj.data[particle, self.obj.lib['y']]),
                                      self.scale['x'](self.obj.data[particle, self.obj.lib['z']]))
        string += 'VERTICES %i %i\n' % (self.obj.particles, self.obj.particles * 2)
        for particle in range(self.obj.particles):
            string += '1 %i\n' % (particle)
        return string

    def property_string(self, prop):
        string = 'SCALARS %s float 1\n' % prop
        string += 'LOOKUP_TABLE default\n'
        for particle in range(self.obj.particles):
            string += self.outlib[prop] % self.scale[prop](self.obj.data[particle, self.obj.lib[prop]])
            string += '\n'
        return string

    def write(self, outputfile, addlib=None, custom_props=None, **kwargs):
        if not addlib is None:
            for lib in addlib:
                self.outlib[lib[0]] = lib[1]
                self.scale[lib[0]] = lambda x: x * lib[2]
        if custom_props is None:
            props = []
            prop = ['r', 'id', 'type', 'aggregate', 'cluster', 'brokenBondFlag', 'percolationPath', 'c_sa', 'c_da', 'c_tot', "relDisplacementLast", "relDisplacementTotal"]
            for p in prop:
                if p in self.obj.lib:
                    props.append(p)
        else:
            props = custom_props

        string = self.header()
        string += self.coordinates()
        string += 'POINT_DATA %i\n' % self.obj.particles
        for prop in props:
            string += self.property_string(prop)
        with open(outputfile, 'w') as f:
            f.write(string)
        boxname = outputfile[:-4] + '_box.vtk'
        with open(boxname, 'w') as f:
            f.write(self.vtk_box())

    def vtk_box(self):
        string = '# vtk DataFile Version 2.0\n'
        string += 'Generated by dda\n'
        string += 'ASCII\n'
        string += 'DATASET RECTILINEAR_GRID\n'
        string += 'DIMENSIONS 2 2 2\n'
        str = ['X', 'Y', 'Z']
        for dim in range(3):
            string += '%s_COORDINATES 2 float\n' % str[dim]
            string += '%e %e\n' % (self.scale['x'](self.obj.box[dim, 0]), self.scale['x'](self.obj.box[dim, 1]))
        return string


import numpy as np


class writedatafile:

    def __init__(self, obj, units='si', style='granular molecular', density=4230, nanodem=True):

        self.obj = obj
        self.units = units
        self._diameters = None
        self.nanodem = nanodem
        if self.units == self.obj.units:
            self.scaleup = 1
        elif self.units == 'si' and self.obj.units == 'nano':
            self.scaleup = 1e-9
        elif self.units == 'micro' and self.obj.units == 'nano':
            self.scaleup = 1e-6

        if style not in ('granular molecular', 'granular bond/gran'):
            raise ValueError('Invalid style, choose from "granular molecular" and "granular bond/gran"')
        self.style = style

        if density is None:
            self.density = float(input('Enter a density for the particles\n> '))
        else:
            self.density = float(density)

    def adjustTypes(self):
        """
        Adjust the types of the particles to get a type for each particle diameter. This is necessary for liggghts bond models
        """
        radii = sorted(list(set(self.obj.data[:, self.obj.lib['r']])))
        for pp in range(self.obj.particles):
            for dia in range(len(radii)):
                if self.obj.data[pp, self.obj.lib['r']] == radii[dia]:
                    self.obj.data[pp, self.obj.lib['type']] = dia + 1
                    break

    @property
    def outlib(self):
        return {'id': '%i',
                'type': '%i',
                'aggregate': '%i',
                'cluster': '%i',
                'x': '%e',
                'y': '%e',
                'z': '%e',
                'r': '%e',
                'diameter': '%e',
                'rho': '%1.3e',
                }

    @property
    def scale(self):
        return {'id': lambda x: x,
                'type': lambda x: x,
                'aggregate': lambda x: x,
                'cluster': lambda x: x,
                'x': lambda x: x * self.scaleup,
                'y': lambda x: x * self.scaleup,
                'z': lambda x: x * self.scaleup,
                'r': lambda x: x * self.scaleup,
                }

    def write(self, outputfile, addlib=None, custom_props=None, **kwargs):
        with open(outputfile, 'w') as f:
            f.write(self.getHeader())
            f.write(self.getGlobalData())
            f.write(self.getBoxData())
            # f.write(self.getMasses())
            f.write(self.getData())

    def getHeader(self):
        return 'LIGGHTS data file (%s, unit system %s)\n\n' % (self.style, self.units)

    def getGlobalData(self):
        string = ''
        if self.style == 'granular bond/gran':
            string += '%i atoms\n' % self.obj.particles
            string += '%i atom types\n' % self.obj.types
            string += '%i bond types\n' % sum(np.arange(self.obj.types))
            string += '0 angle types\n'
            string += '24 extra bond per atom\n\n'
        elif self.style == 'granular molecular':
            string += '%i atoms\n' % self.obj.particles
            string += '%i atom types\n\n' % (self.obj.types + 2)
        return string

    def getBoxData(self):
        string = '%e %e xlo xhi\n' % (self.scale['x'](self.obj.box[0, 0]), self.scale['x'](self.obj.box[0, 1]))
        string += '%e %e ylo yhi\n' % (self.scale['x'](self.obj.box[1, 0]), self.scale['x'](self.obj.box[1, 1]))
        string += '%e %e zlo zhi\n\n' % (self.scale['x'](self.obj.box[2, 0]), self.scale['x'](self.obj.box[2, 1]))
        return string

    def getMasses(self):
        d = np.asarray(list(set(self.diameters))) - 0.15e-9
        masses = 1 / 6 * np.pi * pow(d, 3) * self.density
        string = 'Masses\n\n'
        for r in range(len(masses)):
            string += '%i %e\n' % (r + 1, masses[r])
        string += '\n'
        return string

    def getData(self):
        string = 'Atoms\n\n'
        if self.style in ('granular molecular', 'granular bond/gran'):
            props = ['id', 'type', 'x', 'y', 'z', 'diameter', 'rho', 'aggregate']
            for pp in range(self.obj.particles):
                for prop in props:
                    if prop in ('id', 'type', 'x', 'y', 'z', 'aggregate'):
                        string += self.outlib[prop] % self.scale[prop](self.obj.data[pp, self.obj.lib[prop]])
                    elif prop == 'diameter':
                        string += self.outlib[prop] % self.diameters[pp]
                    elif prop == 'rho':
                        string += self.outlib[prop] % self.density
                    string += ' '
                string += '\n'
        return string

    def getDiameters(self):
        """
        Check if the radius is already adjusted and calculate the diameter
        """
        dmin = 2 * min(self.obj.data[:, self.obj.lib['r']])
        if self.nanodem:
            add = 0.15e-09  # *scaleup
            print('Adding 0.15 nm to diameter')
        else:
            add = 0
        if dmin * 10 - int(dmin * 10) == 0:
            self._diameters = self.scale['r'](2 * (self.obj.data[:, self.obj.lib['r']])) + add

        else:
            self._diameters = self.scale['r'](2 * (self.obj.data[:, self.obj.lib['r']]))  # + 0.15e-09 * self.scaleup

    @property
    def diameters(self):
        if self._diameters is None:
            self.getDiameters()
        return self._diameters

    @diameters.setter
    def diameters(self, array):
        self._diameters = array


class writexyzr:

    def __init__(self, obj, units='nano'):
        self.obj = obj
        self.units = units
        if self.units == self.obj.units:
            self.scaleup = 1
        elif self.units == 'si' and self.obj.units == 'nano':
            self.scaleup = 1e-9
        elif self.units == 'nano' and self.obj.units == 'si':
            self.scaleup = 1e9
        else:
            self.scaleup = 1

    @property
    def outlib(self):
        return {'id': '%i',
                'x': '%e',
                'y': '%e',
                'z': '%e',
                'r': '%e',
                }

    @property
    def scale(self):
        return {'id': lambda x: x,
                'x': lambda x: x * self.scaleup,
                'y': lambda x: x * self.scaleup,
                'z': lambda x: x * self.scaleup,
                'r': lambda x: x * self.scaleup,
                }

    def write(self, outputfile, custom_props=None, **kwargs):
        if custom_props is None:
            props = []
            prop = ['id', 'x', 'y', 'z', 'r']
            for p in prop:
                if p in self.obj.lib:
                    props.append(p)
        else:
            props = []
            for p in custom_props:
                if p in self.outlib:
                    props.append(p)
        string = '{}\n'.format(self.obj.particles)
        string += 'xyzr file in {} format\n'.format(self.units)
        for pp in range(self.obj.particles):
            for prop in props:
                string += self.outlib[prop] % self.scale[prop](self.obj.data[pp, self.obj.lib[prop]])
                string += ' '
            string += '\n'
        with open(outputfile, 'w') as f:
            f.write(string)
