
"""
    This is a module which generates aggregates according to global paramters stated in the main program.

    Usage: create_aggregate(N)

    Source:     Norbert Riefler (University of Bremen)
    Literature: Filippov et al. (2000), Journal of Colloid and Interface Science 229, 261-273
                Zhang et al. (2012), Aerosol Science and Technology, 46: 1065-1078
                Chan and Dahneke (1981), J. Appl. Phys., 52:3106-3110
                W. Hess, H.L. Frisch, R. Klein, Z. Phys. B 64 (1986) 65–67.
    Author:     Valentin Baric (v.baric@iwt.uni-bremen.de)
    Status:     Version 1.0 alpha (02.03.2015)
"""
#############################################
############ Import of Modules ##############
#############################################
import numpy as np
import sys, time, random, math, os
import matplotlib.pyplot as plt
from scipy import stats, spatial
import functions as func
source_file = os.path.dirname(os.path.realpath(__file__))
import datetime
#############################################
######## Generation of an aggregate #########
#############################################
# def create_aggregate(N, rmin, rmax,sigma,TYPE,params):
def create_aggregate(N, rmin, rmax, r_mean, D_f, k_f, epsilon, delta, \
                    deformation, dist, sigma, geom_center, \
                    iteration_limit,debugging, T_gas, TYPE, \
                    rho_1, rho_2, method,dt,mean_r,predef_file,custom_seed,clusters,clustermode,diffusionfile,mixingratio, nAgg, totalAggs):

    if 8 in debugging:
        random.seed(custom_seed)                    # Create a predefined seed (based on the current date)
        np.random.seed(custom_seed)
    else:
        random.seed()                               # Create a new seed (based on the current date)
    ############### Choose Material ############
    if TYPE == 1:
        rho_p = rho_1
    elif TYPE == 2:
        rho_p = rho_2
    restartGlobal = True
    initial = True
    Ninit = np.copy(N)
    clustersizeinit = np.copy(clusters)
    if D_f != 1.0:
        while restartGlobal:
            restartGlobal = False
            ############### Initialising ###############
            j = 0
            loops = 0
            ################ Variables ##################
            r = np.zeros(N)                                 # Radii of the spherules
            x = np.zeros(N)                                 # x-Koordinates of the spherules
            y = np.zeros(N)                                 # y-Koordinates of the spherules
            z = np.zeros(N)                                 # z-Koordinates of the spherules
            cont_i = np.zeros(N)                            # Attempts to set the spherules

                # r_mean = rmin + (rmax-rmin)/2
            # Define the clustersizes
            def get_clustersizes(clustermode, N,clusters):
                """
                Evaluates from the inputfile whether the number of particles per cluster
                or the amount of clusters is defined.
                returns: clusters - <int> the amount of clusters
                         cluster - <array> the particles per cluster (will be returned as empty array here)
                """
                if clustermode == 'CCAN':
                    if clusters > N:
                        clusters = N
                    cluster = np.zeros(clusters,dtype=int)
                elif clustermode == 'CCAS':
                    clusters = int(N/clusters)
                    if clusters == 0:
                        clusters = 1
                    cluster = np.zeros(clusters,dtype=int)
                return clusters,cluster

            # get the amount of <clusters> and the particles per <cluster>
            clusters,cluster = get_clustersizes(clustermode,Ninit,clustersizeinit)
            # This is the amount of clusters of the same size
            index = 0
            cl = 0
            while index < N:
                cluster[cl] += 1
                index += 1
                cl += 1
                if cl == clusters:
                    cl = 0
            types = np.zeros(clusters)

            # Set the type per cluster and shuffle the array to make it random
            if clusters % (1/mixingratio) == 0:
                for index in range(int(len(types)*mixingratio)):
                    types[index] = 1
            else:
                randn = random.randint(0,1)
                if randn == 1:
                    for index in range(int(len(types)*mixingratio) + 1):
                        types[index] = 1
                else:
                    for index in range(int(len(types)*mixingratio)):
                        types[index] = 1


            random.shuffle(types)

            mat_types = np.zeros(N)
            clusterindex = np.zeros(N)
            ###########################################################################
            ###########################################################################
            ###########################################################################
            x_c = []
            y_c = []
            z_c = []
            r_c = []
            index = 0
            for cl in range(clusters):
                xtemp = []
                ytemp = []
                ztemp = []
                rtemp = []
                for cl2 in range(int(cluster[cl])):
                    xtemp.append(0)
                    ytemp.append(0)
                    ztemp.append(0)
                    rtemp.append(0)
                    mat_types[index] = types[cl] + 1
                    clusterindex[index] = cl
                    index += 1
                x_c.append(xtemp)
                y_c.append(ytemp)
                z_c.append(ztemp)
                r_c.append(rtemp)
            particleindex = 0
            clus1 = 0
            while clus1 < clusters:

                PCA = True
                while PCA:
                    PCA = False

                    ############### PP-Distribution ###############
                    if dist == 'lognorm':
                        sigma_dist = np.sqrt(np.log((pow(sigma,2)/pow(rmin,2))+1))
                        mu = np.log(rmin) + pow(sigma_dist,2)/2
                        # use sigma and not sigma_dist here!!!
                        rgen = np.random.lognormal(mu, sigma, cluster[clus1])
                        V_tot = 0
                        A_tot = 0
                        for rr in range(cluster[clus1]):
                            r_c[clus1][rr] = rgen[rr] * pow(10,-9)
                            V_tot += 4/3 * np.pi * pow(r_c[clus1][rr],3)
                            A_tot += 4 * np.pi * pow(r_c[clus1][rr],2)

                    elif dist == 'norm':
                        rgen = np.random.normal(rmin, sigma, cluster[clus1])
                        for rr in range(cluster[clus1]):
                            r_c[clus1][rr] = rgen[rr] * pow(10,-9)
                            V_tot += 4/3 * np.pi * pow(r_c[clus1][rr],3)
                            A_tot += 4 * np.pi * pow(r_c[clus1][rr],2)

                    elif dist == 'none':
                        rmax = rmin
                        for rr in range(cluster[clus1]):
                            r_c[clus1][rr] = rmin
                    elif dist == 'predef':
                        with open(source_file + '/' + predef_file) as f:
                            lines = f.readlines()
                        V_tot = 0
                        A_tot = 0
                        r_prob = np.zeros((len(lines),2))
                        for xxx in range(len(lines)):
                            r_prob[xxx,0] = float(lines[xxx].split()[0])
                            if xxx == 0:
                                r_prob[xxx,1] = float(lines[xxx].split()[1])
                            else:
                                r_prob[xxx,1] = r_prob[xxx-1,1] + float(lines[xxx].split()[1])
                        for xxx in range(cluster[clus1]):
                            r_rand = random.random()
                            xx = 0
                            while r_rand >= r_prob[xx,1]:
                                xx += 1
                            r_c[clus1][xxx] = r_prob[xx,0] * pow(10,-9)
                            V_tot += 4/3 * np.pi * pow(r_c[clus1][xxx],3)
                            A_tot += 4 * np.pi * pow(r_c[clus1][xxx],2)

                    ############### Definition of the first spherule ###############
                    k = 0
                    n = k+1                                             # Number of current spherule
                    # (is always k+1, because k starts at 0)
                    x_c[clus1][k],y_c[clus1][k],z_c[clus1][k] = 0,0,0                              # Place particle at origin
                    ############### Radius of the first spherule ###############
                    if dist == 'random':
                        r_c[clus1][k] = rmin + (rmax - rmin) * random.random()   # Calculate radius within rang
                    max_i = 0                                       # highest attempt to set a particle

                    if cluster[clus1] > 1:                                           # Stop if primary particles are wanted

                        ############### Definition of the second spherule ###############
                        k = 1
                        n = k+1                                         # Particle counter for conveniece
                        if dist == 'random':
                            r_c[clus1][k] = rmin + (rmax - rmin) * random.random()# Random radius deposition
                        rr = epsilon * sum(r_c[clus1][:n])                       # Distance to the next spherule

                        ############### Define a random point on the deposition sphere ###############
                        z_z = rr * (1 - 2 * random.random())            # Random z value
                        theta = np.arccos(z_z/rr)                       # Angle Theta between radius and z - axis
                        phi = random.random() * 2 * np.pi               # Random angle between x- and y axis
                        x_c[clus1][k] = rr * np.cos(phi) * np.sin(theta)         # Resulting x-coordinate
                        y_c[clus1][k] = rr * np.sin(phi) * np.sin(theta)         # Resulting y-coordinate
                        z_c[clus1][k] = z_z                                      # Random z value (see above)
                        rr = 0                                          # reset for next spherule
                        for k in range (2,cluster[clus1]):
                            n = k+1                                     # Particle counter for conveniece
                            # Calculating a new radius if no distribution is used
                            if dist == 'random':
                                r_c[clus1][k] = rmin + (rmax - rmin) * random.random()

                            # --- Calculation of the geometrical center of the current cluster
                            if dist == 'none' or geom_center == 1:
                                # 1) This holds for identical spherules according to [Filippov2000]:
                                x0 = sum(x_c[clus1][:k])/k
                                y0 = sum(y_c[clus1][:k])/k
                                z0 = sum(z_c[clus1][:k])/k
                            else:
                                # 2) These conditions hold for any spherule-radius
                                # Attention must be paid for great size differences. If the first spherule is too big, the center of mass always lies within the first spherule.
                                # => No generation process will be initiated

                                ## The total mass is computed at first (same density, therefore the volume is used)
                                M = sum(np.array(r_c[clus1][:k])**3)
                                x0 = sum(np.array(r_c[clus1][:k])**3 * x_c[clus1][:(n-1)])/M
                                y0 = sum(np.array(r_c[clus1][:k])**3 * y_c[clus1][:(n-1)])/M
                                z0 = sum(np.array(r_c[clus1][:k])**3 * z_c[clus1][:(n-1)])/M

                            # --- Calculation of the position-sphere-radius rr for the next spherule
                            # Used is the Sequental Algorithm (SA) [Filippov2000] Equation [10]:
                            rr = np.sqrt((pow(n,2)*pow(r_mean,2)/(n-1))*((n/k_f)**(2/D_f)) - n*pow(r_mean,2)/(n-1) - n*pow(r_mean,2)*(((n-1)/k_f)**(2/D_f)))
                            # --- Deformate the cluster
                            rr=rr*deformation
                            kk = 0
                            i = 0                                       # reset trial counter
                            # Using the neighborslist to speed up distance search
                            treecoords = []
                            for pp in range(k):
                                treecoords.append( [x_c[clus1][pp],y_c[clus1][pp],z_c[clus1][pp]] )
                            tree = spatial.cKDTree(treecoords)
                            # maximum radius
                            dr = max(r_c[clus1][:k])
                            i = 0
                            done1 = False
                            while not done1:                           # x[k] != 0, if the new position is saved                            # --- Define a random point on the outer sphere
                                z_z = rr * (1 - 2 *random.random())     # Random z value
                                theta = np.arccos(z_z/rr)               # Angle Theta between radius and z - axis
                                phi = random.random() * 2 * np.pi       # Random angle between x- and y axis

                                # --- Random coordinates on Rg Sphere
                                xx = rr*np.cos(phi)*np.sin(theta)
                                yy = rr*np.sin(phi)*np.sin(theta)
                                zz = z_z
                                # --- Shift according to the center of the aggregate
                                xx = xx + x0
                                yy = yy + y0
                                zz = zz + z0

                                c,f = 0,0                               # Counter for try and error
                                # check the amount of particles within reach
                                neighbors = tree.query_ball_point([xx, yy, zz], ( (r_c[clus1][k] + dr) * delta ) )
                                for kk in neighbors:
                                    # --- Calculate the distance
                                    distance = func.distance(x_c[clus1][kk],xx,y_c[clus1][kk],yy,z_c[clus1][kk],zz)
                                    # Check whether the radius rrr is within the section defined by epsilon and delta
                                    if distance > epsilon * (r_c[clus1][k] + r_c[clus1][kk]) and distance < delta * (r_c[clus1][k] + r_c[clus1][kk]):
                                        c += 1
                                    # Check whether the radius rrr is too low
                                    if distance < epsilon * (r_c[clus1][k] + r_c[clus1][kk]):
                                        f += 1

                                if (f == 0 and c>=1):           # Save position of the particle
                                    x_c[clus1][k] = xx
                                    y_c[clus1][k] = yy
                                    z_c[clus1][k] = zz
                                    done1 = True
                                i += 1                          # Trial Counter
                                if i > iteration_limit:         # Too many attempts to set a particle are reached
                                    PCA = True
                                    break                       # Exit "while not done1" loop into for k in range... loop
                                    # The cluster will be rebuild
                            if PCA:
                                break                           # Exit for k in range loop into while PCA loop
                if clus1 == 0 and not PCA:
                    for pp in range(cluster[clus1]):
                        x[particleindex] = x_c[clus1][pp]
                        y[particleindex] = y_c[clus1][pp]
                        z[particleindex] = z_c[clus1][pp]
                        r[particleindex] = r_c[clus1][pp]
                        particleindex += 1
                    clus1 += 1

                elif clus1 != 0 and not PCA:
                    def get_center(x_c,y_c,z_c):
                        x_center = np.mean(x_c)
                        y_center = np.mean(y_c)
                        z_center = np.mean(z_c)
                        return x_center,y_center,z_center

                    def get_shift(x,y,z,x_c,y_c,z_c):
                        xcenter,ycenter,zcenter = get_center(x_c,y_c,z_c)
                        xshift = x - xcenter
                        yshift = y - ycenter
                        zshift = z - zcenter
                        return xshift, yshift, zshift

                    def shiftcluster(x,y,z,x_c,y_c,z_c):
                        xshift, yshift, zshift = get_shift(x,y,z,x_c,y_c,z_c)
                        for pp in range(len(x_c)):
                            x_c[pp] += xshift
                            y_c[pp] += yshift
                            z_c[pp] += zshift

                    def gyrationradius(x,y,z,r,N):
                        M = sum(np.array(r)**3)
                        x0 = sum(np.array(r)**3 * x)/M
                        y0 = sum(np.array(r)**3 * y)/M
                        z0 = sum(np.array(r)**3 * z)/M
                        Rg = 0
                        for kk in range(N):
                            Rg = Rg + pow((x[kk]-x0),2) + pow((y[kk]-y0),2) + pow((z[kk]-z0),2)
                        return np.sqrt(Rg/N)

                    # extract the already set particles
                    x_set = np.zeros(particleindex)
                    y_set = np.zeros(particleindex)
                    z_set = np.zeros(particleindex)
                    r_set = np.zeros(particleindex)

                    for pp in range(particleindex):
                        x_set[pp] = x[pp]
                        y_set[pp] = y[pp]
                        z_set[pp] = z[pp]
                        r_set[pp] = r[pp]
                    N1 = len(x_set)
                    N2 = cluster[clus1]
                    # Shift them so that the geometrical! (see paper) center is in the origin
                    xshift, yshift, zshift = get_shift(0,0,0,x_set,y_set,z_set)
                    for pp in range(len(x_set)):
                        x_set[pp] += xshift
                        y_set[pp] += yshift
                        z_set[pp] += zshift
                    for pp in range(particleindex):
                        x[pp] += xshift
                        y[pp] += yshift
                        z[pp] += zshift
                    # Calculate gyration radii
                    Rg1 = gyrationradius(x_set,y_set,z_set,r_set,N1)
                    Rg2 = gyrationradius(x_c[clus1],y_c[clus1],z_c[clus1],r_c[clus1],N2)
                    broken = False
                    try:
                        rr = np.sqrt((pow(r_mean,2) * pow(N1 + N2,2) / (N1 * N2)) * \
                         pow((N1 + N2)/k_f, (2/D_f)) - \
                         ((N1 + N2) / N2) * pow(Rg1,2) - \
                         ((N1 + N2) / N1) * pow(Rg2,2))
                    except:
                        PCA = True
                        broken = True
                        restartGlobal = True

                    if not broken:
                        # Build a neighbors tree
                        treecoords = []
                        for pp in range(len(x_set)):
                            treecoords.append( [x_set[pp],y_set[pp],z_set[pp]] )
                        tree = spatial.cKDTree(treecoords)
                        # Maximum radius
                        dr = max(r_set)

                        done2 = False
                        i = 0
                        while not done2:
                            # --- Define a random point on the outer sphere
                            z_z = rr * (1 - 2 *random.random())     # Random z value
                            theta = np.arccos(z_z/rr)               # Angle Theta between radius and z - axis
                            phi = random.random() * 2 * np.pi       # Random angle between x- and y axis

                            # --- Random coordinates on Rg Sphere
                            xx = rr*np.cos(phi)*np.sin(theta)
                            yy = rr*np.sin(phi)*np.sin(theta)
                            zz = z_z
                            xshift,yshift,zshift =get_shift(xx,yy,zz,
                                                            x_c[clus1],
                                                            y_c[clus1],
                                                            z_c[clus1])
                            for pp in range(len(x_c[clus1])):
                                x_c[clus1][pp] += xshift
                                y_c[clus1][pp] += yshift
                                z_c[clus1][pp] += zshift

                            c,f = 0,0
                            for pp in range(len(x_c[clus1])):
                                # get all neighbors within radius
                                neighbors = tree.query_ball_point([x_c[clus1][pp],
                                                                   y_c[clus1][pp],
                                                                   z_c[clus1][pp]],
                                                                   ( (r_c[clus1][pp] + dr) * delta ) )
                                # for every neighbor within neighborslist
                                for neigh in neighbors:
                                    # check the distance to the neighbor
                                    distance = func.distance(x_set[neigh],x_c[clus1][pp],
                                                         y_set[neigh],y_c[clus1][pp],
                                                         z_set[neigh],z_c[clus1][pp])
                                    # check if the distance is within threshold (epsilon,delta)
                                    if ((distance > epsilon * (r_set[neigh]+r_c[clus1][pp]))) and \
                                       (distance < delta * (r_set[neigh]+r_c[clus1][pp])):
                                       c += 1
                                    # check if the radius is too low, so there is an overlap
                                    if distance < epsilon * (r_set[neigh]+r_c[clus1][pp]):
                                        f += 1
                            if f == 0 and c >= 1:
                                for pp in range(len(x_c[clus1])):
                                    x[particleindex] = x_c[clus1][pp]
                                    y[particleindex] = y_c[clus1][pp]
                                    z[particleindex] = z_c[clus1][pp]
                                    r[particleindex] = r_c[clus1][pp]
                                    particleindex += 1
                                done2 = True
                                clus1 += 1
                            else:
                                i += 1
                                j += 1
                                if j > 5 * iteration_limit:
                                    restartGlobal = True
                                    break # Out of the while not done loop into while clus1 < cluster loop

                                if i > iteration_limit:
                                    # It is not possible to find a possible spot, generate
                                    # a new cluster and try atain
                                    break # Out of the while not done2 loop into while clus1 < cluster loop
                    if restartGlobal:
                        break
                if restartGlobal:
                    break # Out of for while clus1 < cluster loop into restartGlobal loop
    else:
        # If the fractal dimension is 1.0, then the aggregate equals a straight
        # line. therefore the particles are just put into a line
        loops = 0
        r = np.zeros(N)                                 # Radii of the spherules
        x = np.zeros(N)                                 # x-Koordinates of the spherules
        y = np.zeros(N)                                 # y-Koordinates of the spherules
        z = np.zeros(N)                                 # z-Koordinates of the spherules



        # 2. Define the clustersizes
        # This is the amount of clusters of the same size
        if clusters > N:
            clusters = N
        cluster = np.zeros(clusters,dtype=int)
        index = 0
        cl = 0
        while index < N:
            cluster[cl] += 1
            index += 1
            cl += 1
            if cl == clusters:
                cl = 0
        types = np.zeros(clusters)

        if clusters % 2 == 0:
            for index in range(int(len(types)/2)):
                types[index] = 1
        else:
            randn = random.randint(0,1)
            if randn == 1:
                for index in range(int(len(types)/2 + 1)):
                    types[index] = 1
            else:
                for index in range(int(len(types)/2)):
                    types[index] = 1
        random.shuffle(types)
        mat_types = np.zeros(N)
        clusterindex = np.zeros(N)

        # 1. Primary particle diameters

        if dist == 'lognorm':
            sigma_dist = np.sqrt(np.log((pow(sigma,2)/pow(rmin,2))+1))
            mu = np.log(rmin) + pow(sigma_dist,2)/2
            # use sigma and not sigma_dist here!!!
            r = np.random.lognormal(mu, sigma, N)
            V_tot = 0
            A_tot = 0
            for rr in range(N):
                r[rr] = r[rr] * pow(10,-9)
                V_tot += 4/3 * np.pi * pow(r[rr],3)
                A_tot += 4 * np.pi * pow(r[rr],2)

        elif dist == 'norm':
            r = np.random.normal(rmin, sigma, N)
            for rr in range(N):
                r[rr] = r[rr] * pow(10,-9)
                V_tot += 4/3 * np.pi * pow(r[rr],3)
                A_tot += 4 * np.pi * pow(r[rr],2)

        elif dist == 'none' or 'random':
            rmax = rmin
            for rr in range(cluster[clus1]):
                r[rr] = rmin + (rmax - rmin) * random.random()
        elif dist == 'predef':
            with open(source_file + '/' + predef_file) as f:
                lines = f.readlines()
            V_tot = 0
            A_tot = 0
            r_prob = np.zeros((len(lines),2))
            for xxx in range(len(lines)):
                r_prob[xxx,0] = float(lines[xxx].split()[0])
                if xxx == 0:
                    r_prob[xxx,1] = float(lines[xxx].split()[1])
                else:
                    r_prob[xxx,1] = r_prob[xxx-1,1] + float(lines[xxx].split()[1])
            for xxx in range(N):
                r_rand = random.random()
                xx = 0
                while r_rand >= r_prob[xx,1]:
                    xx += 1
                r[xxx] = r_prob[xx,0] * pow(10,-9)
                V_tot += 4/3 * np.pi * pow(r[xxx],3)
                A_tot += 4 * np.pi * pow(r[xxx],2)


        index = 0
        for cl in range(clusters):
            for cl2 in range(int(cluster[cl])):
                mat_types[index] = types[cl] + 1
                clusterindex[index] = cl
                if index == 0:
                    pass
                else:
                    # Put the particle next to the precedessor
                    x[index] = x[index - 1] + r[index] + r[index -1]
                index += 1

        # rotate the aggregate randomly
        alpha = random.random() * 2 *np.pi
        ca = np.cos(alpha)
        sa = np.sin(alpha)
        C = 1 - ca
        # random vector
        e1 = random.uniform(-10,10)
        e2 = random.uniform(-10,10)
        e3 = random.uniform(-10,10)
        length = np.sqrt( pow(e1, 2) + pow(e2, 2) + pow(e3, 2) )
        e1 /= length
        e2 /= length
        e3 /= length

        # Multiplications
        xs = e1 * sa
        ys = e2 * sa
        zs = e3 * sa
        xC = e1 * C
        yC = e2 * C
        zC = e3 * C
        xyC = e1 * yC
        yzC = e2 * zC
        zxC = e3 * xC

        R_matrix = np.matrix(( (e1*xC+ca,  xyC - zs,  zxC + ys),
                               (xyC + zs,  e2*yC+ca,  yzC - xs),
                               (zxC - ys,  yzC + xs,  e3*zC + ca) ))
        x,y,z = np.array(np.split((R_matrix * [x,y,z]).transpose(),3,axis=1))
        x = np.array([x[0] for x in x])
        y = np.array([y[0] for y in y])
        z = np.array([z[0] for z in z])

    #################################################
    ######## Generation of Aggregate is done ########
    #################Postprocessing##################
    #################################################

    # --- Calculation of the radius of gyration
    n = N
    max_i = 0
    x0 = sum(x[:n])/n
    y0 = sum(y[:n])/n
    z0 = sum(z[:n])/n

    r_g = 0
    for kk in range(n):
        r_g = r_g + pow((x[kk]-x0),2) + pow((y[kk]-y0),2) + pow((z[kk]-z0),2)
    r_g = np.sqrt(r_g/n)

    # --- Calculation of the number density for a cluster
    m = k_f * pow((r_g/r_mean),D_f)
    m = np.round(m)
    # --- Calculation of the smallest radius surrounding the cmplete cluster, centered at the gemetrical center of the cluster:
    rr = np.zeros((n))
    for kk in range(n):
        rr[kk] = np.sqrt( pow((x[kk]-x0),2) + pow((y[kk]-y0),2) + pow((z[kk]-z0),2) )
    rr = max(rr)
    #########################
    #### Diffusion-coeff ####
    #########################
    # ----- Constants ----- #
    k_B = 1.3806488e-23                             # Boltzmann-Constant in J/K
    d_p = r_mean * 2                                # mean diameter in m
    # ----- Calculations ----- #
    rho_gas = 1.20 * 293 / T_gas                    # Density of the gas at temperature T

    # --- kinetic viscosity of air acc. to Daubert 1994 in m2/s
    kin_vis_gas = 1.425e-6 * pow(T_gas,0.5039) / (1 + 1.0883e2 / T_gas) / rho_gas
    # --- mean free path of a gas molecule
    # acc. to Willeke 1976 eq. 18 in m
    lambda_gas = 66e-9 * T_gas / 293 * ((1 + 110.4/293)/(1 + 110.4/T_gas))
    Kn = 2 * lambda_gas / d_p                   # Knudsen-Number
    # --- Calculation of the friction factor
    if method == 1:
        # Chan and Dahneke 1981:
        frict_p = 9.17 * N * kin_vis_gas * (d_p/2) / Kn
        diff_p = k_B * T_gas / frict_p
        # --- Other values which are not further implemented yet
        m_p = N * np.pi / 6 * pow(d_p,3) * rho_p       # Particle mass
        beta = frict_p / m_p
        t_s = pow(r_mean,2) / (2 * diff_p)      # Root Mean Square Displacement (RMSD)
    elif method == 2:
        # Slip correction acc. to Friedlander 2000 eq. 2.21
        Cc = 1 + 2 * lambda_gas / d_p * (1.257 + 0.4 * np.exp(-0.55 * d_p / lambda_gas))
        frict_p = 3 * np.pi * kin_vis_gas * d_p/(Cc)
        diff_p = k_B * T_gas / frict_p
        # --- Other values which are not further implemented yet
        m_p = N * np.pi / 6 * pow(d_p,3) * rho_p       # Particle mass
        beta = frict_p / m_p
        t_s = pow(r_mean,2) / (2 * diff_p)      # Root Mean Square Displacement (RMSD)
    elif method == 3:
        # Calculation method acc to Zhang et al. 2012
        # The hydrodynamic diameter Rh is calculated acc to. Hess et al. (1986)
        diff_data = np.loadtxt(diffusionfile)
        mu,sigma = 0,0
        for lin in range(len(diff_data)):
            if int(diff_data[lin,0]) == N:
                mu = diff_data[lin,1]
                sigma = diff_data[lin,2]
                break
        if mu == 0:
            print('ERROR: %i particles not found in %s' %(N,diffusionfile))
            raise IndexError
        # Return a random value given by a gaussian normal distribution
        # shifted by mu (mean diffusion) and scaled by sigma (std derivation of the diffusion)
        diff_p = sigma * np.random.randn() + mu
        frict_p = k_B * T_gas / diff_p

        m_p = 0
        for pp in range(len(r)):
            m_p += np.pi / 6 * pow(2*r[pp],3) * rho_p
        beta = frict_p / m_p
        t_s = pow(r_mean,2) / (2 * diff_p)




    #########################
    ##### Dimensionless #####
    #########################
    r = r / r_mean
    x = x / r_mean
    y = y / r_mean
    z = z / r_mean

    if dist == 'lognorm' or dist == 'norm':
        r_s = r_mean
    elif dist == 'predef':
        r_s = r_mean
    elif dist == 'none':
        r_s = rmin
        rmin /= rmin
    else:
        r_s = rmin + (rmax-rmin)/2
        rmax = rmax / rmin
        rmin = rmin / rmin

    sys.stdout.write('Aggregate %i/%i complete\n' %(nAgg, totalAggs))
    sys.stdout.flush()
    return N,r,x,y,z,r_g,m,rr,loops,max_i,diff_p,m_p,frict_p,beta,t_s,r_s,Kn,rho_p,mat_types, clusterindex

