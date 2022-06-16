import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def write_liggghts_granular(params, aggs,N_tot, source_file, output_folder, name, BOX, r_ref, max_z):
    INDEX = 1
    outputfile = output_folder + 'data.' + name + '_hybrid_GranMol'

    def read_predef(inputfile):
        """
            Read the file with the predefined diameters
        """
        with open(inputfile) as read:
            lines = read.readlines()
            read.close()
        # Read the predef file and store its content in predef
        predef = np.zeros((len(lines),3))
        for lin in range(len(lines)):
            tmp = lines[lin].split()
            for y in range(len(tmp)):
                predef[lin,y] = tmp[y]
        return predef

    def check_particle_size(predef):
        """
        lists containing the effective diameters and types
        The simulation must not necessarily start with the smallest diameter
        written in the predef file
        """
        effective_dp = []
        effective_type = []
        # check which particle sizes are actually present within the simulation
        for agg in range(len(aggs)):
            for PP in range(aggs[agg].Np):
                i = 0
                while np.around(aggs[agg].r[PP] * r_ref,1) > predef[i,0]:
                    i += 1
                effective_dp.append(predef[i,2]*pow(10,-9))
                effective_type.append(i)
        return effective_dp, effective_type

    def write_existing_types(outputfile,existingType,effectiveList,predef):
        with open(outputfile,'w') as write:
            write.write('#predefindex datafileindex diameter\n')
            for val in range(len(existingType)):
                write.write('%i %i %2.2f\n' %(existingType[val]+1,effectiveList[val,0]+1,predef[int(effectiveList[val,1]),2]))
            write.close()#

    ###############################
    ########### Code states here ##
    ###############################
    # Read the predef file
    predef = read_predef(source_file + '/' + params.pre_dp_file)
    # read the effective diameter and type per particle
    dp_e, type_e = check_particle_size(predef)
    # A list of the types present in the set
    existingType = sorted(list(set(type_e)))
    # change the list type_e so it starts with one and evolves monotoniously
    effectiveList = np.zeros((len(existingType),2))
    for val in range(len(existingType)):
        effectiveList[val,0] = val
        effectiveList[val,1] = existingType[val]
    type_e = [effectiveList[ np.where(effectiveList[:,1] == type_e[val])[0][0],0] for val in range(len(type_e))]
    # Write outputfile
    existing_file = output_folder + 'data.' + name + '_existing_types'
    write_existing_types(existing_file,existingType,effectiveList,predef)

    with open(outputfile, 'w') as write :
        # Head
        write.write('LIGGGHTS data file from DLA Simulation (hybrid granular/molecular)\n\n')
        write.write('%i atoms\n' %N_tot)
        # Atomtypes:
        write.write('%i atom types\n\n' %len(set(type_e)))
        # Box:
        write.write('%1.18e %1.18e xlo xhi\n' %((BOX[0]*r_ref)*pow(10,-9), (BOX[1]*r_ref)*pow(10,-9)))
        write.write('%1.18e %1.18e ylo yhi\n' %((BOX[2]*r_ref)*pow(10,-9), (BOX[3]*r_ref)*pow(10,-9)))
        write.write('0 %1.18e zlo zhi\n\n' %(max_z*r_ref*pow(10,-9)+1*(pow(10,-6))))
        # Masses:
        write.write('Masses\n\n')
        index = 1
        for typ in range(len(effectiveList)):
            # Only the ones which are present in the simulation
            tmp_mass = pow(predef[int(effectiveList[typ,1]),2] * pow(10,-9),3 ) * np.pi * (1/6) * params.rho_1
            write.write('%i %e\n' %(index, tmp_mass))
            index += 1
        #Data for each particle
        write.write('\nAtoms\n\n')
        ind = 0
        for agg in range(len(aggs)):
            for PP in range(aggs[agg].Np):
                write.write('%i %i %1.3e %1.3e %1.12e %1.12e %1.12e %i\n' %(INDEX,
                                                                        type_e[ind]+1,
                                                                        aggs[agg].x[PP]*r_ref*pow(10,-9),
                                                                        aggs[agg].y[PP]*r_ref*pow(10,-9),
                                                                        aggs[agg].z[PP]*r_ref*pow(10,-9),
                                                                        dp_e[ind],
                                                                        aggs[agg].rho,
                                                                        agg+1))
                INDEX += 1
                ind += 1
        write.close()

def PP_histogram(aggs,params,rmin,sigma,output_folder,name):
    ###--- Histogram of primary paricle size
    container_r = []            # Container for radii to calculate size distribution
    for k in range(len(aggs)):
        for size in range(len(aggs[k].r)):
            container_r.append(aggs[k].r[size]*params.rmin)

    fig = plt.figure(figsize=(10,5))
    ax = fig.add_subplot(111)
    count, bins, ignored = ax.hist(container_r, 100, align='mid')

    if params.r_dist == 'lognorm':
        sigma_dist = np.sqrt(np.log((pow(sigma,2)/pow(rmin,2))+1))
        mu = np.log(rmin) - pow(sigma_dist,2)/2
        x = np.linspace(min(bins), max(bins), 10000)
        pdf = (np.exp(-(np.log(x) - mu)**2 / (2 * sigma**2)) / (x * sigma * np.sqrt(2 * np.pi)))
        plt.plot(x, pdf, linewidth=2, color='r')
        plt.axis('tight')
        #ax.set_xscale('log')
        ax.set_xlim(1,50)
    elif params.r_dist == 'norm':
        plt.plot(bins, 1/(sigma * np.sqrt(2 * np.pi)) * \
                       np.exp( - (bins - rmin)**2 / (2 * sigma**2) ),
                       linewidth=2, color='r')
        upper_limit = max(container_r)*1.1
        ax.set_xlim(0,upper_limit)
    else:
        upper_limit = max(container_r)*1.1
        ax.set_xlim(0,upper_limit)
    ax.set_ylim(0,max(count)*1.1)
    ax.set_xlabel('Radius / nm')
    ax.set_ylabel('counts')
    plt.tight_layout()
    plt.savefig(output_folder + name + '_PP_histogram' + params.img_format)
    plt.close()

def N_histogram(N_PP_per_Agg,output_folder,name,params):
    fig = plt.figure(figsize=(10,5))
    ax = fig.add_subplot(111)
    N = max(N_PP_per_Agg) - min(N_PP_per_Agg) +1
    count, bins, ignored = ax.hist(N_PP_per_Agg,N, align='mid')
    ax.set_ylim(0,max(count)*1.1)
    ax.set_xlim(0,max(N_PP_per_Agg)*1.1)
    ax.set_xlabel('Primary Particles per Aggregate / -')
    ax.set_ylabel('counts')
    plt.tight_layout()
    plt.savefig(output_folder + name + '_N_PP_histogram' + params.img_format)
    plt.close()
