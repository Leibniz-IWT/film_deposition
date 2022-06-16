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
#############################################
######## Generation of an aggregate #########
#############################################
# def create_aggregate(N, rmin, rmax,sigma,TYPE,params):
def create_aggregate(N, rmin, rmax, r_mean, D_f, k_f, epsilon, delta, \
                    deformation, dist, sigma, geom_center, \
                    iteration_limit,debugging, T_gas, TYPE, \
                    rho_1, rho_2, method,dt,mean_r,predef_file,custom_seed):
    ############### Choose Material ############
    if TYPE == 1:
        rho_p = rho_1
    elif TYPE == 2:
        rho_p = rho_2
    mat_type = [TYPE for pp in range(N)]
    clusterindex = np.zeros(N)
    ############### Initialising ###############
    loops = 0                                           # Attempts to generate an aggregate
    restart = 1                                         # Restart trigger will be set if program is stuck
    while restart == 1:                                 # Restart loop
        ################ Variables ##################
        if 8 in debugging:
            random.seed(custom_seed)                    # Create a predefined seed (based on the current date)
            np.random.seed(custom_seed)
        else:
            random.seed()                               # Create a new seed (based on the current date)
        r = np.zeros((N))                               # Radii of the spherules
        x = np.zeros((N))                               # x-Koordinates of the spherules
        y = np.zeros((N))                               # y-Koordinates of the spherules
        z = np.zeros((N))                               # z-Koordinates of the spherules
        cont_i = np.zeros((N))                              # Attempts to set the spherules

        ############### PP-Distribution ###############
        if dist == 'lognorm':
            sigma_dist = np.sqrt(np.log((pow(sigma,2)/pow(rmin,2))+1))
            mu = np.log(rmin) + pow(sigma_dist,2)/2
            # use sigma and not sigma_dist here!!!
            r = np.random.lognormal(mu, sigma, N)
            V_tot = 0
            A_tot = 0
            for rr in range(len(r)):
                r[rr] = r[rr] * pow(10,-9)
                V_tot += 4/3 * np.pi * pow(r[rr],3)
                A_tot += 4 * np.pi * pow(r[rr],2)
            # ----- Calculation of mean diameter
            # There are several aequivalent diameters:
            #   0: This is the value used for the distrobution (normal: arth/lognormal: geom)
            #   1: The arithmetic mean of the PP radii distribution
            #   2: The geometric mean of the PP radii distribution
            #   3: The harmonic mean of the PP radii (no real relevance)
            #   4: The sauter diameter of the PP radii distribution
            # if mean_r == 0:
            #     r_mean = rmin * pow(10,-9)
            # elif mean_r == 1:
            #     r_mean = np.mean(r)
            # elif mean_r == 2:
            #     r_mean = stats.mstats.gmean(r)
            # elif mean_r == 3:
            #     r_mean = stats.mstats.hmean(r) # mean radius in m
            # elif mean_r == 4:
            #     r_mean = 3* V_tot / A_tot
        elif dist == 'norm':
            r = np.random.normal(rmin, sigma, N)
            for rr in range(len(r)):
                r[rr] = r[rr] * pow(10,-9)
                V_tot += 4/3 * np.pi * pow(r[rr],3)
                A_tot += 4 * np.pi * pow(r[rr],2)
            # if mean_r == 0:
            #     r_mean = rmin * pow(10,-9)
            # elif mean_r == 1:
            #     r_mean = np.mean(r)
            # elif mean_r == 2:
            #     r_mean = stats.mstats.gmean(r)
            # elif mean_r == 3:
            #     r_mean = stats.mstats.hmean(r)
            # elif mean_r == 4:
            #     r_mean = 3* V_tot / A_tot
        elif dist == 'none':
            #r_mean = rmin
            rmax = rmin
        elif dist == 'predef':
            read = open(source_file + '/' + predef_file)
            lines = read.readlines()
            read.close()
            V_tot = 0
            A_tot = 0
            r_prob = np.zeros((len(lines),2))
            for xxx in range(len(lines)):
                r_prob[xxx,0] = float(lines[xxx].split()[0])
                if xxx == 0:
                    r_prob[xxx,1] = float(lines[xxx].split()[1])
                else:
                    r_prob[xxx,1] = r_prob[xxx-1,1] + float(lines[xxx].split()[1])
            for xxx in range(len(r)):
                r_rand = random.random()
                xx = 0
                while r_rand >= r_prob[xx,1]:
                    xx += 1
                r[xxx] = r_prob[xx,0] * pow(10,-9)
                V_tot += 4/3 * np.pi * pow(r[xxx],3)
                A_tot += 4 * np.pi * pow(r[xxx],2)
            # if mean_r == 0:
            #     r_mean = rmin
            # elif mean_r == 1:
            #     r_mean = np.mean(r)
            # elif mean_r == 2:
            #     r_mean = stats.mstats.gmean(r)
            # elif mean_r == 3:
            #     r_mean = stats.mstats.hmean(r)
            # elif mean_r == 4:
            #     r_mean = 3* V_tot / A_tot
        else:
            pass
            # r_mean = rmin + (rmax-rmin)/2
        if D_f != 1.0:

            ############### Definition of the first spherule ###############
            restart = 0
            k = 0
            n = k+1                                             # Number of current spherule
            # (is always k+1, because k starts at 0)
            x[k],y[k],z[k] = 0,0,0                              # Place particle at origin
            ############### Radius of the first spherule ###############
            if dist == 'none' or dist == 'random':
                r[k] = rmin + (rmax - rmin) * random.random()   # Calculate radius within rang
            max_i = 0                                       # highest attempt to set a particle

            if N > 1:                                           # Stop if primary particles are wanted
                ############### Definition of the second spherule ###############
                k = 1
                n = k+1                                         # Particle counter for conveniece
                if dist == 'none' or dist == 'random':
                    r[k] = rmin + (rmax - rmin) * random.random()# Random radius deposition
                rr = epsilon * sum(r[:n])                       # Distance to the next spherule

                ############### Define a random point on the deposition sphere ###############
                z_z = rr * (1 - 2 * random.random())            # Random z value
                theta = np.arccos(z_z/rr)                       # Angle Theta between radius and z - axis
                phi = random.random() * 2 * np.pi               # Random angle between x- and y axis
                x[k] = rr * np.cos(phi) * np.sin(theta)         # Resulting x-coordinate
                y[k] = rr * np.sin(phi) * np.sin(theta)         # Resulting y-coordinate
                z[k] = z_z                                      # Random z value (see above)

                rr = 0                                          # reset for next spherule
                i = 0                                           # trial counter
                for k in range (2,N):
                    if i > max_i:
                        max_i = i                               # Counter for maximal attempts
                    n = k+1                                     # Particle counter for conveniece
                    # Calculating a new radius if no distribution is used
                    if dist == 'none' or dist == 'random':
                        r[k] = rmin + (rmax - rmin) * random.random()

                    # --- Calculation of the geometrical center of the current cluster
                    if dist == 'none' or geom_center == 1:
                        # 1) This holds for identical spherules according to [Filippov2000]:
                        x0 = sum(x[:(n-1)])/(n-1)
                        y0 = sum(y[:(n-1)])/(n-1)
                        z0 = sum(z[:(n-1)])/(n-1)
                    else:
                        # 2) These conditions hold for any spherule-radius
                        # Attention must be paid for great size differences. If the first spherule is too big, the center of mass always lies within the first spherule.
                        # => No generation process will be initiated

                        ## The total mass is computed at first (same density, therefore the volume is used)
                        M = sum(r[:(n-1)]**3)
                        x0 = sum(r[:(n-1)]**3 * x[:(n-1)])/M
                        y0 = sum(r[:(n-1)]**3 * y[:(n-1)])/M
                        z0 = sum(r[:(n-1)]**3 * z[:(n-1)])/M

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
                        treecoords.append( [x[pp],y[pp],z[pp]] )
                    tree = spatial.cKDTree(treecoords)
                    # maximum radius
                    dr = max(r[:k])
                    while x[k] == 0:                           # x[k] != 0, if the new position is saved
                        # --- Define a random point on the outer sphere
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
                        rrr = np.zeros((k))
                        # check the amount of particles within reach
                        neighbors = tree.query_ball_point([xx, yy, zz], ( (r[k] + dr) * delta ) )
                        for kk in neighbors:
                            # --- Calculate the distance
                            rrr[kk] = func.distance(x[kk],xx,y[kk],yy,z[kk],zz)
                            # Check whether the radius rrr is within the section defined by epsilon and delta
                            if ((rrr[kk] > epsilon * (r[k]+r[kk]))) and (rrr[kk] < delta*(r[k]+r[kk])):
                                c=c+1
                            # Check whether the radius rrr is too low
                            if rrr[kk] < epsilon*(r[k]+r[kk]):
                                f=f+1

                        if (f == 0 and c>=1):           # Save position of the particle
                            x[k] = xx
                            y[k] = yy
                            z[k] = zz
                            cont_i[k] = i
                        i += 1                          # Trial Counter
                        if i > iteration_limit:         # Too many attempts to set a particle are reached
                            restart = 1                 # Set trigger to restart the generation process
                            break                       # Exit "while [x[k]=0]" loop

                    if restart == 1:                    # Call if restart escape was called
                        loops += 1                      # How many restarts couter
                        break                           # exit for loop
        else:
            restart = 0
            for pp in range(N):
                ############### Radius of the spherule ###############
                if dist == 'none' or dist == 'random':
                    r[pp] = rmin + (rmax - rmin) * random.random()   # Calculate radius within rang
                if pp == 0:
                    pass
                else:
                    x[pp] = x[pp - 1] + r[pp] + r[pp - 1]

    if D_f == 1.0:
        #-------Trig.-functions-and-axis-------
        alpha = random.random() * 2 *np.pi
        ca = np.cos(alpha)
        sa = np.sin(alpha)
        C = 1 - ca
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

        R_matrix = np.matrix(((e1*xC + ca, xyC - zs, zxC + ys),(xyC + zs, e2*yC + ca, yzC - xs),(zxC - ys, yzC + xs, e3*zC + ca) ))
        inputmatrix = np.zeros((N, 3))
        outputmatrix = np.zeros((N, 3))

        #------calculate-coordinates-----------
        for i in range(N):
            for j in range(0,3):
                outputmatrix[i, j] += R_matrix[j, 0] * x[i]
                outputmatrix[i, j] += R_matrix[j, 1] * y[i]
                outputmatrix[i, j] += R_matrix[j, 2] * z[i]
        for pp in range(N):
            x[pp] = outputmatrix[pp,0]
            y[pp] = outputmatrix[pp,1]
            z[pp] = outputmatrix[pp,2]

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
    m = round(m)

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
    elif method == 2:
        # Slip correction acc. to Friedlander 2000 eq. 2.21
        Cc = 1 + 2 * lambda_gas / d_p * (1.257 + 0.4 * np.exp(-0.55 * d_p / lambda_gas))
        frict_p = 3 * np.pi * kin_vis_gas * d_p/(Cc)
    elif method == 3:
        # Calculation method acc to Zhang et al. 2012
        # The hydrodynamic diameter Rh is calculated acc to. Hess et al. (1986)
        Rh = (2*(D_f - 1))/(D_f) * rr
        PA = pow(N/0.505,1.03) * pow(r_mean,2)*np.pi
        frict_p = (6*np.pi*kin_vis_gas*Rh)/(1+(lambda_gas*np.pi*Rh/PA)*(1.257+0.4*np.exp(-1.1*PA/(lambda_gas*np.pi*Rh))))

    diff_p = k_B * T_gas / frict_p
    # --- Other values which are not further implemented yet
    m_p = N * np.pi / 6 * pow(d_p,3) * rho_p       # Particle mass
    beta = frict_p / m_p
    t_s = pow(r_mean,2) / (2 * diff_p)      # Root Mean Square Displacement (RMSD)
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
    return N,r,x,y,z,r_g,m,rr,loops,max_i,diff_p,m_p,frict_p,beta,t_s,r_s,Kn,rho_p,mat_type, clusterindex

