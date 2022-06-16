import numpy as np
from scipy import stats
import sys, random

def calc_r(params,source_file):
    rmin = params.rmin * pow(10,-9)
    if params.r_dist == 'none':
        rmin = params.rmin * pow(10,-9)
        rmax = rmin
        sigma = 0
        r_mean = rmin
    elif params.r_dist == 'random':
        rmin = params.rmin * pow(10,-9)
        rmax = params.rmax * pow(10,-9)
        sigma = 0
        r_mean = rmin + (rmax-rmin)/2
    elif params.r_dist == 'lognorm':
        rmin = params.rmin
        rmax = params.rmin
        sigma = params.r_sigma
        r_mean = rmin * pow(10,-9)
    elif params.r_dist == 'norm':
        rmin = params.rmin
        rmax = params.rmin
        sigma = params.r_sigma
        r_mean = rmin * pow(10,-9)
    elif params.r_dist == 'predef':
        rmin = params.rmin * pow(10,-9)
        rmax = params.rmin * pow(10,-9)
        sigma = 0
        read = open(source_file + '/' + params.pre_dp_file)
        lines = read.readlines()
        read.close()
        rList = np.zeros(len(lines))
        VTot = 0
        ATot = 0
        for xxx in range(len(lines)):
            rList[xxx] = float(lines[xxx].split()[0]) *pow(10,-9)
            VTot += 4/3 * np.pi * pow(rList[-1],3)
            ATot += 4 * np.pi * pow(rList[-1],2)
        if params.mean_r == 0:
            r_mean = rmin
        elif params.mean_r == 1:
            r_mean = np.mean(rList)
        elif params.mean_r == 2:
            r_mean = stats.mstats.gmean(rList)
        elif params.mean_r == 3:
            r_mean = stats.mstats.hmean(rList)
        elif params.mean_r == 4:
            r_mean = 3* VTot / ATot
    else:
        raise NotImplementedError("""
            r_dist in inputfile is not valid
            """)
    return rmin,rmax,sigma,r_mean
    # normal and log-normal distributions are further handled within the creation module

def calc_N(params, source_file, N_tot):
    # ----- Amount of Particles per aggregate --- #
    N_PP_per_Agg = []           # Container for PP per Agg before the generation
    N_PP = 0                    # Number of current particles
    if params.N_dist == 'none' or params.N_dist == 'random':
        while N_PP < N_tot:
            """ Prepare the amount of PP per Agg. This has to be done before the generation to use parallel computing"""
            if params.N_dist == 'none':
                N = int(params.N_min)
            elif params.N_dist == 'random':
                N = int(params.N_min + (params.N_max - params.N_min) * random.random())
            N_PP_per_Agg.append(N)
            N_PP += N
    elif params.N_dist == 'norm' or params.N_dist == 'lognorm':
        if params.N_dist == 'lognorm':
            sigma_N_dist =np.sqrt(np.log((pow(params.N_sigma,2)/pow(params.N_min,2))+1))
            N_mu = np.log(params.N_min) - pow(sigma_N_dist,2)/2
            N_distr = np.random.lognormal(N_mu,sigma_N_dist,int((N_tot / params.N_min) * 2))
        elif params.N_dist == 'norm':
            N_distr = np.random.normal(params.N_min, params.N_sigma, (N_tot / params.N_min) * 2)
        i = 0
        while N_PP < N_tot:
            N_PP_per_Agg.append(int(N_distr[i]))
            N_PP += int(N_distr[i])
            i += 1
    elif params.N_dist == 'predef':
        read = open(source_file + '/' + params.pre_N_file)
        lines = read.readlines()
        read.close()
        N_prob = np.zeros((len(lines),2))
        for xxx in range(len(lines)):
            N_prob[xxx,0] = int(lines[xxx].split()[0])
            if xxx == 0:
                N_prob[xxx,1] = float(lines[xxx].split()[1])
            else:
                N_prob[xxx,1] = N_prob[xxx-1,1] + float(lines[xxx].split()[1])
        while N_PP < N_tot:
            N_rand = random.random()
            xx = 0
            while N_rand >= N_prob[xx,1]:
                xx += 1
            N_PP_per_Agg.append(int(N_prob[xx,0]))
            N_PP += int(N_prob[xx,0])
    else:
        sys.exit('Check your inputfile, N_dist is invalid')
    # ----- If the maximal number of particles is exceeded
    if N_PP > N_tot and params.N_dist != 'predef':
        N_PP_per_Agg[-1] -= N_PP - N_tot
        N_PP = N_tot
    elif N_PP > N_tot and params.N_dist == 'predef':
        # Remove the last aggregate in case a predef distribution with
        # distince amount of particles is wanted
        N_PP -= N_PP_per_Agg[-1]
        N_tot = N_PP
        N_PP_per_Agg.pop()
    return N_PP, N_tot, N_PP_per_Agg
