import sys
import os
from writelayer import writeVTK
from writelayer import writetrj
from writelayer import writedatafile
from writelayer import writexyzr
from scipy import spatial
import numpy as np


class layer:

    def __init__(self):

        self._particles = None
        self._data = None
        self._box = None
        self._lines = None
        self._lib = None
        self._units = None
        self._filetype = None
        self._types = None
        self._coordination = None
        self._kdtree = None
        self._meancoordination = None
        self._sumcoordination = None
        self._rmax = None
        self._xtree = None
        self._ytree = None
        self._ztree = None
        self._rtree = None
        self._aggtree = None
        self._idtree = None
        self.file = None
        self._rbmax = None
        self._addlib = None
        self._nAggregates = None
        self._totalparticlevolume = None

    def write(self, outputfile, units='nano', custom_props=None, ** kwargs):
        if self.getFormat(outputfile) == 'vtk':
            output = writeVTK(self, units)
        elif self.getFormat(outputfile) == 'trj':
            output = writetrj(self, units)
        elif self.getFormat(outputfile) == 'data':
            output = writedatafile(self, units)
        elif self.getFormat(outputfile) == 'xyzr':
            output = writexyzr(self, units)
        output.write(outputfile, self._addlib, custom_props=custom_props, **kwargs)

    def getFormat(self, filename):
        if filename.endswith('vtk'):
            return 'vtk'
        elif filename.endswith('trj') or 'dump' in filename:
            return 'trj'
        elif '.data' in filename or \
             'data.' in filename:
            return 'data'
        elif filename.endswith('xyzr'):
            return 'xyzr'
        else:
            raise ValueError('Unknown outputfile format')

    @property
    def particles(self):
        """
            :return: The amount of particles within the box
        """
        return self._particles

    @particles.setter
    def particles(self, amount):
        if type(amount) == int:
            self._particles = amount
        else:
            raise ValueError('Only integers can be used for amount of particles')

    @property
    def lib(self):
        return self._lib

    @lib.setter
    def lib(self, library):
        if not type(library) == dict:
            raise ValueError('Only libraries are allowed for lib')
        self._lib = library

    def addlib(self, name, stringOut, scale=1):
        """
            Add a custom value to the outputvalues. This requires adding a column to the array of the data
        """
        if name not in self.lib:
            if self._addlib is None:
                self._addlib = []
            self._addlib.append([name, stringOut, scale])
            tmpData = np.zeros((np.shape(self.data)[0], np.shape(self.data)[1] + 1))
            tmpData[:, :-1] = self.data
            self.data = tmpData
            self.lib[name] = np.shape(self.data)[1] - 1
            print('Added property {} to dataset'.format(name))

    @property
    def nAggregates(self):
        """
            Return the number of aggregates in the film
        """
        if self._nAggregates is None:
            self._nAggregates = len(list(set(self.aggregates)))
        return self._nAggregates

    @property
    def aggregates(self):
        """
            This returns a list of the corresponding aggregates of the particles
        """
        return self.data[:, self.lib["aggregate"]]

    @property
    def cluster(self):
        """
            This returns a list of the corresponding cluster of the particles
        """
        return self.data[:, self.lib["cluster"]]

    @property
    def x(self):
        """
            Returns the x values of the layer
        """
        return self.data[:, self.lib["x"]]

    @x.setter
    def x(self, array):
        """
            Returns the x values of the layer
        """
        self._data[:, self.lib["x"]] = array

    @property
    def y(self):
        """
            Returns the y values of the layer
        """
        return self.data[:, self.lib["y"]]

    @y.setter
    def y(self, array):
        """
            Returns the y values of the layer
        """
        self._data[:, self.lib["y"]] = array

    @property
    def z(self):
        """
            Returns the z values of the layer
        """
        return self.data[:, self.lib["z"]]

    @z.setter
    def z(self, array):
        """
            Returns the z values of the layer
        """
        self._data[:, self.lib["z"]] = array

    @property
    def type(self):
        """
            Returns the z values of the layer
        """
        return self.data[:, self.lib["type"]]

    @property
    def r(self):
        """
            Returns the r values of the layer
        """
        return self.data[:, self.lib["r"]]

    @r.setter
    def r(self, array):
        """
            Returns the r values of the layer
        """
        self._data[:, self.lib["r"]] = array

    @property
    def rmax(self):
        """
        Return the maximum radius of self.r
        """
        if self._rmax is None:
            self._rmax = max(self.r)
        return self._rmax

    @rmax.setter
    def rmax(self, value):
        """
        Return the maximum radius of self.r
        """
        self._rmax = value

    @property
    def data(self):
        if self._data is None:
            self.getData()
            # self._adjustIntData()
        return self._data

    @data.setter
    def data(self, matrix):
        self._data = matrix

    def getData(self):
        pass

    def sortData(self, property):
        self._data = self.data[self.data[:, self.lib[property]].argsort(), :]

    @property
    def units(self):
        return self._units

    @units.setter
    def units(self, units):
        if units not in ('nano', 'si'):
            raise ValueError('Unit System not supported')
        self._units = units

    @property
    def box(self):
        if self._box is None:
            self.getBox()
        return self._box

    @box.setter
    def box(self, size):
        if np.shape(size) != (3, 2):
            raise ValueError('Box size does not fit!')
        self._box = size

    def getBox(self):
        pass

    @property
    def filetype(self):
        return self._filetype

    @filetype.setter
    def filetype(self, type):
        self._filetype = type

    @property
    def types(self):
        if self._types is None:
            self._types = len(set((self.data[:, self.lib['type']])))
        return self._types

    @types.setter
    def types(self, value):
        self._types = value

    def _adjustIntData(self):
        if min(self.data[:, self.lib['id']]) != 0:
            for pp in range(self.particles):
                for prop in ('id', 'type', 'aggregate'):
                    self.data[pp, self.lib[prop]] -= 1
                if self.filetype == 'dda':
                    self.data[pp, self.lib['cluster']] -= 1

    def distance(self, p1, p2, squared=True):
        """
        Function that returns the distance of the two particles p1 and p2
        :p1: Index of Particle 1
        :p2: Index of Particle 2
        :squared: If True, the distance will be the squared distance. Omitting the squareroot is faster
        """
        dist = pow(self.data[p1, self.lib['x']] - self.data[p2, self.lib['x']], 2) + \
            pow(self.data[p1, self.lib['y']] - self.data[p2, self.lib['y']], 2) + \
            pow(self.data[p1, self.lib['z']] - self.data[p2, self.lib['z']], 2)
        if not squared:
            dist = np.sqrt(dist)
        return dist

    def distanceTree(self, p1, p2, squared=True):
        """
        Function that returns the distance of the two particles p1 and p2
        :p1: Index of Particle 1
        :p2: Index of Particle 2
        :squared: If True, the distance will be the squared distance. Omitting the squareroot is faster
        """
        dist = pow(self.data[p1, self.lib['x']] - self.xtree[p2], 2) + \
            pow(self.data[p1, self.lib['y']] - self.ytree[p2], 2) + \
            pow(self.data[p1, self.lib['z']] - self.ztree[p2], 2)
        if not squared:
            dist = np.sqrt(dist)
        return dist

    def neighborslist(self, force=False, periodic=True, verbose=False):
        """
        This functions builds a neighborslist out of self.data
        For this the scipy.spatial.cKDTree function is used
        :return: KDTree with the neighborslist
        """
        if self._kdtree is None or force:
            x = self.x
            y = self.y
            z = self.z
            r = self.r
            agg = self.data[:, self.lib['aggregate']]
            idstore = self.data[:, self.lib['id']]
            addwidth = 2 * self.rmax
            totalwidth = (self.box[0, 1] - self.box[0, 0])
            totalInefficient = False
            count = 0
            # print(addwidth, totalwidth, self.box, self.x[0])
            if periodic:
                if totalInefficient:
                    for particle in range(self.particles):
                        x = np.append(x, self.x[particle] + totalwidth)
                        y = np.append(y, self.y[particle])
                        z = np.append(z, self.z[particle])
                        r = np.append(r, self.r[particle])
                        agg = np.append(agg, self.data[particle, self.lib['aggregate']])
                        idstore = np.append(idstore, self.data[particle, self.lib['id']])

                        x = np.append(x, self.x[particle] - totalwidth)
                        y = np.append(y, self.y[particle])
                        z = np.append(z, self.z[particle])
                        r = np.append(r, self.r[particle])
                        agg = np.append(agg, self.data[particle, self.lib['aggregate']])
                        idstore = np.append(idstore, self.data[particle, self.lib['id']])

                        x = np.append(x, self.x[particle])
                        y = np.append(y, self.y[particle] + totalwidth)
                        z = np.append(z, self.z[particle])
                        r = np.append(r, self.r[particle])
                        agg = np.append(agg, self.data[particle, self.lib['aggregate']])
                        idstore = np.append(idstore, self.data[particle, self.lib['id']])

                        x = np.append(x, self.x[particle])
                        y = np.append(y, self.y[particle] - totalwidth)
                        z = np.append(z, self.z[particle])
                        r = np.append(r, self.r[particle])
                        agg = np.append(agg, self.data[particle, self.lib['aggregate']])
                        idstore = np.append(idstore, self.data[particle, self.lib['id']])

                        x = np.append(x, self.x[particle] + totalwidth)
                        y = np.append(y, self.y[particle] + totalwidth)
                        z = np.append(z, self.z[particle])
                        r = np.append(r, self.r[particle])
                        agg = np.append(agg, self.data[particle, self.lib['aggregate']])
                        idstore = np.append(idstore, self.data[particle, self.lib['id']])

                        x = np.append(x, self.x[particle] + totalwidth)
                        y = np.append(y, self.y[particle] - totalwidth)
                        z = np.append(z, self.z[particle])
                        r = np.append(r, self.r[particle])
                        agg = np.append(agg, self.data[particle, self.lib['aggregate']])
                        idstore = np.append(idstore, self.data[particle, self.lib['id']])

                        x = np.append(x, self.x[particle] - totalwidth)
                        y = np.append(y, self.y[particle] + totalwidth)
                        z = np.append(z, self.z[particle])
                        r = np.append(r, self.r[particle])
                        agg = np.append(agg, self.data[particle, self.lib['aggregate']])
                        idstore = np.append(idstore, self.data[particle, self.lib['id']])

                        x = np.append(x, self.x[particle] - totalwidth)
                        y = np.append(y, self.y[particle] - totalwidth)
                        z = np.append(z, self.z[particle])
                        r = np.append(r, self.r[particle])
                        agg = np.append(agg, self.data[particle, self.lib['aggregate']])
                        idstore = np.append(idstore, self.data[particle, self.lib['id']])
                else:
                    for particle in range(self.particles):
                        if self.x[particle] < self.box[0, 0] + addwidth:
                            x = np.append(x, self.x[particle] + totalwidth)
                            y = np.append(y, self.y[particle])
                            z = np.append(z, self.z[particle])
                            r = np.append(r, self.r[particle])
                            agg = np.append(agg, self.data[particle, self.lib['aggregate']])
                            idstore = np.append(idstore, self.data[particle, self.lib['id']])
                        elif self.x[particle] > self.box[0, 1] - addwidth:
                            x = np.append(x, self.x[particle] - totalwidth)
                            y = np.append(y, self.y[particle])
                            z = np.append(z, self.z[particle])
                            r = np.append(r, self.r[particle])
                            agg = np.append(agg, self.data[particle, self.lib['aggregate']])
                            idstore = np.append(idstore, self.data[particle, self.lib['id']])
                        if self.y[particle] < self.box[1, 0] + addwidth:
                            x = np.append(x, self.x[particle])
                            y = np.append(y, self.y[particle] + totalwidth)
                            z = np.append(z, self.z[particle])
                            r = np.append(r, self.r[particle])
                            agg = np.append(agg, self.data[particle, self.lib['aggregate']])
                            idstore = np.append(idstore, self.data[particle, self.lib['id']])
                        elif self.y[particle] > self.box[1, 1] - addwidth:
                            x = np.append(x, self.x[particle])
                            y = np.append(y, self.y[particle] - totalwidth)
                            z = np.append(z, self.z[particle])
                            r = np.append(r, self.r[particle])
                            agg = np.append(agg, self.data[particle, self.lib['aggregate']])
                            idstore = np.append(idstore, self.data[particle, self.lib['id']])
                        if self.x[particle] < self.box[0, 0] + addwidth and \
                           self.y[particle] < self.box[1, 0] + addwidth:
                            x = np.append(x, self.x[particle] + totalwidth)
                            y = np.append(y, self.y[particle] + totalwidth)
                            z = np.append(z, self.z[particle])
                            r = np.append(r, self.r[particle])
                            agg = np.append(agg, self.data[particle, self.lib['aggregate']])
                            idstore = np.append(idstore, self.data[particle, self.lib['id']])
                        elif self.x[particle] < self.box[0, 0] + addwidth and \
                                self.y[particle] > self.box[1, 1] - addwidth:
                            x = np.append(x, self.x[particle] + totalwidth)
                            y = np.append(y, self.y[particle] - totalwidth)
                            z = np.append(z, self.z[particle])
                            r = np.append(r, self.r[particle])
                            agg = np.append(agg, self.data[particle, self.lib['aggregate']])
                            idstore = np.append(idstore, self.data[particle, self.lib['id']])
                        elif self.x[particle] > self.box[0, 1] - addwidth and \
                                self.y[particle] < self.box[1, 0] + addwidth:
                            z = np.append(z, self.z[particle])
                            r = np.append(r, self.r[particle])
                            agg = np.append(agg, self.data[particle, self.lib['aggregate']])
                            x = np.append(x, self.x[particle] - totalwidth)
                            y = np.append(y, self.y[particle] + totalwidth)
                            idstore = np.append(idstore, self.data[particle, self.lib['id']])
                        elif self.x[particle] > self.box[0, 1] - addwidth and \
                                self.y[particle] > self.box[1, 1] - addwidth:
                            x = np.append(x, self.x[particle] - totalwidth)
                            y = np.append(y, self.y[particle] - totalwidth)
                            z = np.append(z, self.z[particle])
                            r = np.append(r, self.r[particle])
                            agg = np.append(agg, self.data[particle, self.lib['aggregate']])
                            idstore = np.append(idstore, self.data[particle, self.lib['id']])
            data = [x, y, z]
            if verbose:
                print('Neighbortree Length:', len(x))
            self._xtree = x
            self._ytree = y
            self._ztree = z
            self._rtree = r
            self._aggtree = agg
            self._idtree = idstore
            self._kdtree = spatial.cKDTree(np.asarray(data).T.tolist())
        return self._kdtree

    @property
    def xtree(self):
        return self._xtree

    @property
    def rtree(self):
        return self._rtree

    @property
    def ytree(self):
        return self._ytree

    @property
    def ztree(self):
        return self._ztree

    @property
    def aggtree(self):
        return self._aggtree

    @property
    def idtree(self):
        return self._idtree

    def calcCoordination(self, method=1, delta=1.01, threshold=0.5, force=False, periodic=True, verbose=False):
        """
            This calculates the coordination of the particles within the films
            PRECONDITION: This uses the neighborslist (self.neighborslist) and the self.rmax functions

            The neighborslist is used to only include the particles within the range of particle 1's radius + rmax
            This drastically speeds up the calculation

            Methods are:
                1 = (r1 + r2) * delta
                2 = r1 + r2 + (r1 + r2) / 2
                3 = r1 + r2 + threshold
            for the radius which the distance of the two particles are compared to

            Returns nothing

        """
        # Only if it is not calculated yet
        self._rbmax = np.zeros(4)
        self._coordination = np.zeros((self.particles, 3))
        tree = self.neighborslist(force=True, periodic=periodic, verbose=verbose)
        for particle in range(len(self.data)):
            self.data[particle, self.lib['c_sa']] = 0
            self.data[particle, self.lib['c_da']] = 0
            self.data[particle, self.lib['c_tot']] = 0
            coordinate = [self.data[particle, self.lib['x']],
                          self.data[particle, self.lib['y']],
                          self.data[particle, self.lib['z']]]
            neighbors = tree.query_ball_point(coordinate, (self.data[particle, self.lib['r']] + 2 * self.rmax))
            for particle2 in neighbors:
                # The same particle is included in the output of the query_ball_point function
                if particle2 != particle:
                    distance = self.distanceTree(particle, particle2, squared=True)
                    if method == 1:
                        distance_compare = pow((self.data[particle, self.lib['r']] + self.rtree[particle2]) * delta, 2)
                    elif method == 2:
                        distance_compare = pow((self.data[particle, self.lib['r']] + self.rtree[particle2])
                                               + (self.data[particle, self.lib['r']] + self.rtree[particle2]) / 2, 2)
                    elif method == 3:
                        distance_compare = pow((self.data[particle, self.lib['r']] + self.rtree[particle2] + threshold), 2)
                    if distance < distance_compare:
                        if self.data[particle, self.lib['aggregate']] == self.aggtree[particle2]:
                            self._coordination[particle, 0] += 1
                            self.data[particle, self.lib['c_sa']] += 1
                            rb = 0.48 * min([self.rtree[particle], self.rtree[particle2]])
                            if rb > self._rbmax[0]:
                                self._rbmax[0] = rb
                                self._rbmax[1] = np.sqrt(distance)
                                self._rbmax[2] = self.rtree[particle]
                                self._rbmax[3] = self.rtree[particle2]
                        else:
                            self._coordination[particle, 1] += 1
                            self.data[particle, self.lib['c_da']] += 1
                        self.data[particle, self.lib['c_tot']] += 1
        self._coordination[:, 2] = self._coordination[:, 0] + self._coordination[:, 1]

    @property
    def rbmax(self):
        """
            Maximum bond radius
        """
        if self._rbmax is None:
            self.calcCoordination(delta=1.1, force=True)
        return self._rbmax

    @property
    def coordination(self):
        if self._coordination is None:
            self.calcCoordination()
        return self._coordination

    @coordination.setter
    def coordination(self, array):
        self._coordination = array

    @property
    def meancoordination(self):
        if self._meancoordination is None:
            self._meancoordination = np.zeros(3)
            self._meancoordination[0] = np.mean(self.data[:, self.lib['c_sa']])
            self._meancoordination[1] = np.mean(self.data[:, self.lib['c_da']])
            self._meancoordination[2] = self._meancoordination[0] + self._meancoordination[1]
        return self._meancoordination

    @meancoordination.setter
    def meancoordination(self, vector):
        self._meancoordination = vector

    @property
    def sumcoordination(self):
        """
            Calculate the sum of contacts
        """
        if self._sumcoordination is None:
            self._sumcoordination = np.zeros(3)
            self._sumcoordination[0] = np.sum(self.data[:, self.lib['c_sa']])
            self._sumcoordination[1] = np.sum(self.data[:, self.lib['c_da']])
            self._sumcoordination[2] = self._sumcoordination[0] + self._sumcoordination[1]
        return self._sumcoordination

    def printOpenFoamBoxString(self):
        string = "\n"
        string += "\t({:+3.2f} {:+3.2f} {:+3.2f})\n".format(self.box[0, 0], self.box[1, 0], self.box[2, 0])
        string += "\t({:+3.2f} {:+3.2f} {:+3.2f})\n".format(self.box[0, 1], self.box[1, 0], self.box[2, 0])
        string += "\t({:+3.2f} {:+3.2f} {:+3.2f})\n".format(self.box[0, 1], self.box[1, 1], self.box[2, 0])
        string += "\t({:+3.2f} {:+3.2f} {:+3.2f})\n".format(self.box[0, 0], self.box[1, 1], self.box[2, 0])
        string += "\t({:+3.2f} {:+3.2f} {:+3.2f})\n".format(self.box[0, 0], self.box[1, 0], self.box[2, 1])
        string += "\t({:+3.2f} {:+3.2f} {:+3.2f})\n".format(self.box[0, 1], self.box[1, 0], self.box[2, 1])
        string += "\t({:+3.2f} {:+3.2f} {:+3.2f})\n".format(self.box[0, 1], self.box[1, 1], self.box[2, 1])
        string += "\t({:+3.2f} {:+3.2f} {:+3.2f})\n".format(self.box[0, 0], self.box[1, 1], self.box[2, 1])
        print(string)

    def getMinZR(self):
        zr = self.data[:, self.lib['z']] - self.data[:, self.lib['r']]
        zrmin = min(zr)
        return zrmin

    def getMinDim(self, dim):
        return min(self.data[:, self.lib[dim]])

    def getMaxDim(self, dim):
        return max(self.data[:, self.lib[dim]])

    def shiftToZ(self, z):
        zrmin = self.getMinZR()
        shift = z - zrmin
        for pp in range(self.particles):
            self.data[pp, self.lib['z']] += shift

    def adjustbox(self):
        """
            Adjust the box dimensions to match the structure
        """
        print('Adjusting Box')
        self.box[2, 0] = min(self.data[:, self.lib['z']] - self.data[:, self.lib['r']])
        self.box[2, 1] = max(self.data[:, self.lib['z']] + self.data[:, self.lib['r']])

    def getZHistLimits(self, cutoff=10000, scaling=None, stepsize=1e-9):
        """
            Make a histogram over the z-axis and cut off the border region

        """
        self.sortData('z')
        """
        Get the limits for the binning of the histogram.
        """
        minx = min(self.data[:, self.lib['z']])
        maxx = max(self.data[:, self.lib['z']])
        minx2 = min(self.data[cutoff:-cutoff, self.lib['z']])
        maxx2 = max(self.data[cutoff:-cutoff, self.lib['z']])

        bins = np.arange(minx, maxx, stepsize)
        bins2 = np.arange(minx2, maxx2, stepsize)
        hist2, bin_edge2 = np.histogram(self.data[cutoff:-cutoff, self.lib['z']], bins=bins2)
        hist, bin_edge = np.histogram(self.data[:, self.lib['z']], bins=bins)

        if scaling is None:
            scale = 1.00
        else:
            scale = scaling

        mean = np.mean(hist2) * scale
        start = 0
        stop = len(self.data) - 1
        for i in range(len(hist)):
            if hist[i] < mean:
                start = i
            else:
                break
        for i in reversed(range(0, len(hist))):
            if hist[i] > mean:
                stop = i
                break

        return bin_edge[start], bin_edge[stop]

    def removeBoundaryLayers(self, cutoff=10000, scaling=None, stepsize=1e-9):
        """
            Remove the border particles to remove boundary effects

        """
        limits = self.getZHistLimits(cutoff=cutoff, scaling=scaling, stepsize=stepsize)
        indices = []
        for pp in range(self.particles):
            if self.data[pp, self.lib['z']] > limits[1] or \
                    self.data[pp, self.lib['z']] < limits[0]:
                indices.append(pp)
        self._data = np.delete(self.data, indices, 0)
        self._particles = len(self._data)

    def porosity(self, lower=0.1, upper=0.9):
        """
        Calculate the porosity of the film between lower and upper limit
        """
        def scale(value):
            """
            Scale to SI units
            """
            if value < 1e-3:
                return value
            else:
                return value * 1e-9

        top = max(self.z + self.r)
        bot = min(self.z - self.r)
        height = top - bot
        botlim = bot + lower * height
        toplim = top - (1 - upper) * height
        particleVolume = 0
        for pp in range(self.particles):
            if self.z[pp] < botlim or \
                    self.z[pp] > toplim:
                pass
            else:
                particleVolume += 4 / 3 * np.pi * pow(scale(self.r[pp]), 3)
        boxVolume = scale(toplim - botlim) * scale(self.box[0, 1] - self.box[0, 0]) * scale(self.box[1, 1] - self.box[1, 0])
        return 1 - particleVolume / boxVolume

    @property
    def totalparticlevolume(self):
        """
            Calculate the total volume of the particles
        """
        if self._totalparticlevolume is None:
            if self.r[0] > 1e-3:
                scale = 1e-9
            else:
                scale = 1
            self._totalparticlevolume = sum(4 / 3 * np.pi * pow(self.r * scale, 3))
        return self._totalparticlevolume

    def adjustRadiusLIGGGHTS(self):
        """
            Remove the extra particle radius that comes from LIGGGHTS for the water layer
        """
        self.r -= 0.075
