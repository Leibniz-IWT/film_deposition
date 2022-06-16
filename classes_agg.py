################################################
############### Classes ########################
################################################
"""
Definition of the classes Aggregate and Particle which contain the data
for each aggregate and particle, respectively
"""
import random
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt

class Particle:
    """ The relevant data if needed for a single PP """
    def __init__(self,index,agg,agg_pp,r,TYPE):
        self.idx = index            # global index
        self.agg = agg              # agg index
        self.TYPE = TYPE            # Type of material
        self.pp = agg_pp            # PP index within aggregate
        self.r = r                  # radius of the PP
        self.idxtop = self.idx      # Closest neighbor in z-direction (pos)
        self.idxbot = self.idx      # Closest neighbor in z-direction (neg)
        self.x = 0                  # x-position
        self.y = 0                  # y-position
        self.z = 0                  # z-position

        # Percolation parameters
        self.adj = np.zeros(100,dtype=int)
        self.adj_num = 0
        self.wall_contact = 0
        self.path_chain = 0
        self.predec = 0

        # Coordination parameter
        self.c_0 = np.zeros(3,dtype=int)

    def show(self):
        print('Index ', self.idx)
        print('Top Neighbor ', self.idxtop)
        print('Bottom Neighbor ', self.idxbot)
        print('x ', self.x)
        print('y ', self.y)
        print('z ', self.z)
class Aggregate:
    """The data from the create_aggregate function is stored here.
    This data will be contained and never changed"""
    def __init__(self,values,r_dist):

        # PP and Radii
        self.Np = values[0]         # PP in the aggregate
        self.r = values[1]          # radii
        self.dist = r_dist          # PP size distribution method
        # The global index of each PP
        self.index = np.zeros((self.Np,3))
        # Coordinates
        self.x = values[2]          # x-coordinates
        self.x_s = np.zeros((self.Np,1))        # original coordinates
        self.y = values[3]          # y-coordinates
        self.y_s = np.zeros((self.Np,1))        # original coordinates
        self.z = values[4]          # z-coordinates
        self.z_s = np.zeros((self.Np,1))        # original coordinates

        # Statistics
        self.r_g = values[5]        # gyration-radius
        self.nd = values[6]         # Number-density
        self.r_small = values[7]    # smalles sphere around agg
        self.loops = values[8] + 1  # Loops it took to generate Aggregate
        self.max_i = values[9] + 1  # Maxmial attemps to set a particle
        self.D = values[10]         # Diffusion coefficient
        self.m_p = values[11]       # Mass
        self.f = values[12]         # friction factor
        self.beta = values[13]      # Relaxation time
        self.t_s = values[14]       # Scaling time (RMSD)
        self.r_s = values[15]       # scaling length
        self.Kn = values[16]        # Knudsen Number
        self.rho = values[17]       # Density of material that was used
        self.TYPE = values[18]      # Material type index
        self.cluster = values[19]   # Index of the Cluster
        self.c_0 = np.zeros((self.Np,3))
        # Closest neigbours
        # For each particle the closest neighbor in z direction is stored:
        # neigh_high = [Ind_agg_above, Ind_PP_above]
        # neigh_low = [Ind_agg_underneath, Ind_PP_underneath]
        #self.neigh_high = np.zeros((self.Np,3))
        #self.neigh_low = np.zeros((self.Np,3))

        # Booleans which will be changed during the code
        self.cross = False          # Has the aggregate crossed the walls
        self.agg_hit = False        # Boolean for the deposition
        self.init = False           # Boolean for the initiation

        # Coordinates which are needed for the position after a hit
        self.hit = np.zeros((self.Np,1))
    ##########################################################################
    ##########################################################################
    # Function to print all the data
    def show(self):
        print('\n#--- Reading particle data ---#')
        print('Primary Particles:\t\t%s' %self.Np)
        print('Radii (nm):\t\t\t', end='')
        for ii in range(len(self.r)):
            print('%1.3e '%self.r[ii], end = ' ')
        print('\nx-coordinates (m):\t\t', end='')
        for ii in range(len(self.x)):
            print('%1.3e '%self.x[ii], end = ' ')
        print('\ny-coordinates (m):\t\t', end='')
        for ii in range(len(self.y)):
            print('%1.3e '%self.y[ii], end = ' ')
        print('\nz-coordinates (m):\t\t', end='')
        for ii in range(len(self.z)):
            print('%1.3e '%self.z[ii], end = ' ')
        print('\nGyration-Radius (m):\t\t%1.3e' %self.r_g)
        print('Smalles radius (m):\t\t%1.3e' %self.r_small)
        print('Number density:\t\t\t%i' %self.nd)
        print('Mass (kg):\t\t\t%2.3e' %self.m_p)
        print('Diffusion-coefficient (m²/s):\t%2.3e' %self.D)
        print('Friction-Coefficient:\t\t%2.3e' %self.f)
        print('Knudsen-Number:\t\t\t%2.3e' %self.Kn)
        print('Loops to create Particle:\t%i' %self.loops)
        print('Maximum attempts to set PP:\t%i' %self.max_i)
        print('Density (kg/m³):\t\t%i' %self.rho)
        print('#--- End ---#')

    ##########################################################################
    ##########################################################################
    def center(self):
        """ Returns the coordinates of the center """
        x0 = sum(self.x)/len(self.x)
        y0 = sum(self.y)/len(self.y)
        z0 = sum(self.z)/len(self.z)
        return x0,y0,z0

    ##########################################################################
    ##########################################################################
    def print_center(self):
        """ Prints the coordinates of the center """
        x0,y0,z0 = self.center()
        print('\n#--- Center of the aggregate ---#')
        print('x0:\t%1.3e' %x0)
        print('y0:\t%1.3e' %y0)
        print('z0:\t%1.3e' %z0)

    ##########################################################################
    ##########################################################################
    def shift_to(self,x,y,z):
        """ Center an aggregate at a distinct coordinate"""
        x0,y0,z0 = self.center()
        for ii in range(len(self.x)):
            self.x[ii] += x - x0
            self.y[ii] += y - y0
            self.z[ii] += z - z0

    ##########################################################################
    ##########################################################################
    def shift_by(self,x,y,z):
        """ Shift aggregate by x,y,z-range """
        for ii in range(len(self.x)):
            self.x[ii] += x
            self.y[ii] += y
            self.z[ii] += z

    ##########################################################################
    ##########################################################################
    def rotate_rnd(self):
        """ Rotate the aggregate randomly around its center """
        a = 5
        # Get current coodinates
        x0,y0,z0 = self.center()
        r_g0 = 0
        for kk in range(len(self.x)):
            r_g0 = r_g0 + pow((self.x[kk]-x0),2) + pow((self.y[kk]-y0),2) + pow((self.z[kk]-z0),2)
        r_g0 = np.sqrt(r_g0/len(self.x))

        # --- Calculation of the number density for a cluster
        m0 = 1.8 * pow((r_g0/a),1.8)
        m0 = round(m0)

        # --- Calculation of the smallest radius surrounding the cmplete cluster, centered at the gemetrical center of the cluster:
        rr0 = np.zeros((len(self.x)))
        for kk in range(len(self.x)):
            rr0[kk] = np.sqrt( pow((self.x[kk]-x0),2) + pow((self.y[kk]-y0),2) + pow((self.z[kk]-z0),2) )
        rr0 = max(rr0)
        # Prepare rotation matrix
        #-------Trig.-functions-and-axis-------
        # Random angle and axis
        alpha = random.random() * 2 * np.pi
        axis = [random.randint(1,9), random.randint(1,9), random.randint(1,9)]
        len_axis = np.sqrt( axis[0]**2 + axis[1]**2 + axis[2]**2 )
        axis_first = axis
        for jj in range(len(axis)):
            axis[jj] = axis[jj] / len_axis
        # shortenings
        ca = np.cos(alpha)
        sa = np.sin(alpha)
        C = 1 - ca
        x, y, z = (axis)
        xs = x*sa
        ys = y*sa
        zs = z*sa
        xC = x*C
        yC = y*C
        zC = z*C
        xyC = x*yC
        yzC = y*zC
        zxC = z*xC
        # The rotation matrix (wikipedia)
        R_matrix = np.matrix(((x*xC + ca, xyC - zs, zxC + ys),(xyC + zs, y*yC + ca, yzC - xs),(zxC - ys, yzC + xs, z*zC + ca) ))
        # The aggregate is rotated around 0,0,0, therefore it has to be shifted there
        self.shift_to(0,0,0)

        # Generation of a matrix (easier to calculate)
        inputmatrix = np.zeros((len(self.x), 3))
        outputmatrix = np.zeros((len(self.x), 3))
        for ii in range(len(self.x)):
            inputmatrix[ii,0] = self.x[ii]
            inputmatrix[ii,1] = self.y[ii]
            inputmatrix[ii,2] = self.z[ii]
        #inputmatrix = inputmatrix.T
        #print('\ninput.T\n')
        #print(inputmatrix)
        #------calculate-coordinates-----------
        for i in range(len(inputmatrix)):
           for j in range(0,3):
               for n in range(0,3):
                   outputmatrix[i, j] += R_matrix[j, n] * inputmatrix[i,n]
        #outputmatrix = R_matrix * inputmatrix
        #outputmatrix = outputmatrix.T
        #---- Change coordinates in aggregate
        for ii in range(len(self.x)):
            self.x[ii] = outputmatrix[ii,0]
            self.y[ii] = outputmatrix[ii,1]
            self.z[ii] = outputmatrix[ii,2]
        self.shift_to(x0,y0,z0)
        r_g = 0
        for kk in range(len(self.x)):
            r_g = r_g + pow((self.x[kk]-x0),2) + pow((self.y[kk]-y0),2) + pow((self.z[kk]-z0),2)
        r_g = np.sqrt(r_g/len(self.x))

        # --- Calculation of the number density for a cluster
        m = 1.8 * pow((r_g/a),1.8)
        m = round(m)

        # --- Calculation of the smallest radius surrounding the cmplete cluster, centered at the gemetrical center of the cluster:
        rr = np.zeros((len(self.x)))
        for kk in range(len(self.x)):
            rr[kk] = np.sqrt( pow((self.x[kk]-x0),2) + pow((self.y[kk]-y0),2) + pow((self.z[kk]-z0),2) )
        rr = max(rr)
        # no return needed, the data is stored within the class

    ##########################################################################
    ##########################################################################
    def plot(self):
        """ Function to plot the aggregate"""
        fig = plt.figure()
        self.shift_to(0,0,0)
        r = self.r
        x = self.x
        y = self.y
        z = self.z
        # ----- Histogram ----- #
        ax = fig.add_subplot(1,2,1)
        count, bins, ignored = ax.hist(r, 100, normed=True, align='mid')
        if self.dist == 'lognorm':
            ax.set_xscale('log')
            ax.set_xlim(1e-1,1e2)
        else:
            upper_limit = max(r)*1.1
            ax.set_xlim(0,upper_limit)
        ax.set_xlabel('Radius / nm')
        ax.set_ylabel('Fraction / %')
        # ----- Spheres ----- #
        ax = fig.add_subplot(1,2,2, projection='3d')
        ax.set_aspect('equal')
        u = np.linspace(0, 2 * np.pi ,100)
        v = np.linspace(0, np.pi, 100)
        xp,yp,zp = [],[],[]
        N = len(r)
        for i in range(N):
            xp = r[i] * np.outer(np.cos(u), np.sin(v)) + x[i]
            yp = r[i] * np.outer(np.sin(u), np.sin(v)) + y[i]
            zp = r[i] * np.outer(np.ones(np.size(u)), np.cos(v)) + z[i]
            ax.plot_surface(xp, yp, zp, rstride=10, cstride=10, color='k', linewidth=0)
        max_range = np.array([(max(x)-min(x)),(max(y)-min(y)),(max(z)-min(z))]).max()
        LIMIT = max_range
        ax.set_xlim(-LIMIT,LIMIT)
        ax.set_ylim(-LIMIT,LIMIT)
        ax.set_zlim(-LIMIT,LIMIT)
        ax.set_xlabel('x')
        ax.set_ylabel('y')
        ax.set_zlabel('z')
        plt.show()

    ##########################################################################
    ##########################################################################
    def save_plot(self,filename):
        """ Function to save a plot the aggregate"""
        fig = plt.figure()
        self.shift_to(0,0,0)
        r = self.r
        x = self.x
        y = self.y
        z = self.z
        # ----- Histogram ----- #
        ax = fig.add_subplot(1,2,1)
        count, bins, ignored = ax.hist(r, 100, normed=True, align='mid')
        if self.dist == 'lognorm':
            ax.set_xscale('log')
            ax.set_xlim(1e-10,1e-7)
        else:
            upper_limit = max(r)*1.1
            ax.set_xlim(0,upper_limit)
        ax.set_xlabel('Radius / m')
        ax.set_ylabel('Fraction / %')
        # ----- Spheres ----- #
        ax = fig.add_subplot(1,2,2, projection='3d')
        ax.set_aspect('equal')
        u = np.linspace(0, 2 * np.pi ,100)
        v = np.linspace(0, np.pi, 100)
        xp,yp,zp = [],[],[]
        N = len(r)
        for i in range(N):
            xp = r[i] * np.outer(np.cos(u), np.sin(v)) + x[i]
            yp = r[i] * np.outer(np.sin(u), np.sin(v)) + y[i]
            zp = r[i] * np.outer(np.ones(np.size(u)), np.cos(v)) + z[i]
            ax.plot_surface(xp, yp, zp, rstride=10, cstride=10, color='k', linewidth=0)
        max_range = np.array([(max(x)-min(x)),(max(y)-min(y)),(max(z)-min(z))]).max()
        LIMIT = max_range * 1.1
        ax.set_xlim(-LIMIT,LIMIT)
        ax.set_ylim(-LIMIT,LIMIT)
        ax.set_zlim(-LIMIT,LIMIT)
        ax.set_xlabel('x')
        ax.set_ylabel('y')
        ax.set_zlabel('z')
        plt.savefig(filename)

    ##########################################################################
    ##########################################################################
    # ---- Each Particle will be indexed like:
    # [glob_index, agg_index, PP_internal_index]
    # Furthermore the neighbors are set to be itself
    def init_index(self,agg_index,glob_part):
        for PP in range(self.Np):
            self.index[PP] = [glob_part, agg_index, PP]
            self.neigh_high[PP] = self.index[PP]
            self.neigh_low[PP] = self.index[PP]
            glob_part += 1
        return glob_part
    ##########################################################################
    def store_coords(self):
        for PP in range(self.Np):
            self.x_s[PP] = self.x[PP]
            self.y_s[PP] = self.y[PP]
            self.z_s[PP] = self.z[PP]
    ##########################################################################
    def load_coords(self):
        for PP in range(self.Np):
            self.x[PP] = self.x_s[PP]
            self.y[PP] = self.y_s[PP]
            self.z[PP] = self.z_s[PP]

    def checkWall(self, box):

        for pp in range(self.Np):
            if  self.x[pp] < box[0]:
                self.x[pp] -= 2 * box[0]
            elif self.x[pp] > box[1]:
                self.x[pp] -= 2 * box[1]

            if self.y[pp] < box[2]:
                self.y[pp] -= 2 * box[2]
            elif self.y[pp] > box[3]:
                self.y[pp] -= 2 * box[3]
