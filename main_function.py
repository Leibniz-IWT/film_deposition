# --- Python integrated modules
import sys, random, os, datetime, warnings, socket, argparse
import multiprocessing as mp
# --- Third party modules
import numpy as np
from scipy import spatial, stats
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import shutil, fcntl, termios, struct
# --- Own modules
import functions as f
from create_aggregate import create_aggregate as SA
from create_aggregate_CCA import create_aggregate as CCA
import debugging as debug
# This is the aggregate class which stores all relevant data
from classes_agg import Aggregate as Aggregate
from classes_agg import Particle as Particle
import write_output
import terminal_print
import statistics
warnings.filterwarnings('error')
if sys.stdout.isatty():
    Redirection = 0
    COLS = struct.unpack('hh',  fcntl.ioctl(sys.stdout, termios.TIOCGWINSZ, '1234'))[1] -1
else:
    COLS = 60
    Redirection = 1

def init_params(params_file, source_file):
    cwd = os.getcwd()
    if cwd != source_file:
        f_params = cwd + '/'+ params_file
        target = source_file + '/' + params_file
        try:
            shutil.copy2(f_params, target)
        except FileNotFoundError:
            sys.exit('''

                      Error 01: File %s not found. Check directory!

                      '''%params_file )
    if(params_file.endswith(".py")):
        params_file = params_file[:-3]


    exec('import ' + params_file +' as params')


def output_trj(params):
    output_folder = params.output_folder
    PID = os.getpid()  # Process-ID of the current simulation
    now = datetime.datetime.now().replace(microsecond=0)
    # Potential Output of the trajectory
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)
    outputfile = output_folder + params.experiment + params.traj_format
    write_lammpstrj = open(outputfile, 'w')

    if 1 in params.output:
        if not os.path.exists(output_folder):
            os.makedirs(output_folder)
        logfile = output_folder + 'log.' + params.experiment
        log = open(logfile, 'w')
        log.write('Experiment:\t\t%s (PID: %i)\n' % (params.experiment, PID))
        log.write('\t\t\t')
        log.write(now.strftime("%Y-%m-%d %H:%M"))
        log.write('\n')
        log.write('Particles:\t\t%i\n' % params.N_tot)
        log.write('\n##### Initiating Particle Generation')
        log.flush()
    return write_lammpstrj

def generate_film(params, source_file, num_cores, jobserver):
    PID = os.getpid()  # Process-ID of the current simulation
    now = datetime.datetime.now().replace(microsecond=0)
    aggs = []  # Container for all aggregates
    AGG_TYPE = []  # Container for a random distribution of aggregates types

    # Get the amount of particles per aggregate
    rmin, rmax, sigma, r_mean = statistics.calc_r(params, source_file)
    N_PP, N_tot, N_PP_per_Agg = statistics.calc_N(params, source_file, params.N_tot)
    # ----- Material selection
    for PP in range(len(N_PP_per_Agg)):
        if 4 not in params.modules:
            # No distribution
            AGG_TYPE.append(1)
        else:
            if params.Material_selection == 1:
                # Switch for each aggregate
                if PP % 2 == 0:
                    AGG_TYPE.append(1)
                elif PP % 2 != 0:
                    AGG_TYPE.append(2)
            elif params.Material_selection == 2:
                # random distribution
                AGG_TYPE.append(random.randint(1, 2))

    # Print the terminal output with all the data
    terminal_print.loghead(params, params.experiment, PID, num_cores, now, N_PP, N_PP_per_Agg, rmin, rmax)

    #
    if params.DLA_Mode == 'PCA':
        results = [jobserver.apply_async(SA, args=( \
            N_PP_per_Agg[PP], rmin, rmax, r_mean, params.D_f, params.k_f,
            params.epsilon, params.delta, params.deformation,
            params.r_dist, sigma, params.geom_center, params.iteration_limit,
            params.debugging, params.T, AGG_TYPE[PP], params.rho_1,
            params.rho_2, params.Fric_method, params.dt , params.mean_r,
            params.pre_dp_file, params.custom_seed_aggs, params.diffusionfile)) \
                   for PP in range(len(N_PP_per_Agg))]
    if 'CCA' in params.DLA_Mode:
        results = [jobserver.apply_async(CCA, args=( \
            N_PP_per_Agg[PP], rmin, rmax, r_mean, params.D_f, params.k_f,
            params.epsilon, params.delta, params.deformation,
            params.r_dist, sigma, params.geom_center, params.iteration_limit,
            params.debugging, params.T, AGG_TYPE[PP], params.rho_1,
            params.rho_2, params.Fric_method, params.dt , params.mean_r,
            params.pre_dp_file, params.custom_seed_aggs, params.Cluster, params.DLA_Mode, params.diffusionfile,
            params.mixingratio, PP, len(N_PP_per_Agg))) \
                   for PP in range(len(N_PP_per_Agg))]

    results = [p.get() for p in results]
    for res in range(len(results)):
        aggs.append(Aggregate(results[res], params.r_dist))

    return results, aggs, N_PP_per_Agg

def global_pp_indexing(aggs):
    pt = []  # Container which includes the particles
    # ----- Global indexing of each PP ----- #
    # In each instance of Aggregate the global index will be stored
    idx = 0
    # Get the maximal radius of the particles to govern size effects in the
    # neigborlist distance search (set_new_pos)
    maximal_radius = 0
    for agg in range(len(aggs)):
        for PP in range(aggs[agg].Np):
            pt.append(Particle(idx,agg,PP,aggs[agg].r[PP],aggs[agg].TYPE[PP]))
            idx += 1
            if aggs[agg].r[PP] > maximal_radius:
                maximal_radius = aggs[agg].r[PP]
    generation_complete = datetime.datetime.now().replace(microsecond=0)
    return pt, maximal_radius


def particle_transport(aggs, pt, params, maximal_radius, write_lammpstrj):
    # Shift all aggregates initially to a position where they cannot be hit
    for agg in range(len(aggs)):
        aggs[agg].shift_to(0, 0, -2e-7)
    f.ref_pt_pos(pt, aggs)
    count_steps = 1  # The overall steps that were done

    # ----- Box Size Calculation ----- #
    L_x_1 = - params.Box_size / 2
    L_x_2 = params.Box_size / 2
    L_y_1 = - params.Box_size / 2
    L_y_2 = params.Box_size / 2
    Height_above = params.Height_above
    BOX = [L_x_1, L_x_2, L_y_1, L_y_2]
    max_z = 0
    if 1 not in params.debugging:
        # Code 6 in debugging means that there is no movement, just generation of aggs
        print('\n##### Particle Transport calculation initiated')
        sys.stdout.flush()
        max_z = 1 * params.rmin / params.rmin
        np_box = 0  # Particles in the box
        count_t = 0
        count_agg = 0  # aggregates in the box
        count_PP = 0
        av_Pe = []
        av_Vel = []
        iteration_step = 1
        building_tree_time = []
        for k in range(len(aggs)):
            if 9 not in params.debugging and Redirection == 0:
                # Progress bar
                sys.stdout.write(f.progress(k, len(aggs)) + '\r')
                sys.stdout.flush()
            if 1 in params.output and count_steps % params.log_step == 0:
                now_log = datetime.datetime.now().replace(microsecond=0).strftime("%x  %X")
                log.write('%4i (%4i/%4i)\t%s\n' % ((count_steps, k + 1, len(aggs), now_log)))
                log.flush()
            ###################################
            ######## Precalculations ##########
            ###################################
            '''
            The data which is the same for all iterations for the same
            aggregate are calulated once and used until the aggregate
            is deposited
            '''
            ## ----- Parameters of current aggregate
            # stored in short-letter variables for convenience
            r = aggs[k].r  # Radii
            x = aggs[k].x  # Current x-coordinates
            y = aggs[k].y  # Current y-coordinates
            z = aggs[k].z  # Current z-coordinates
            D = aggs[k].D  # Diffusion coefficient of current aggregate
            m_p = aggs[k].m_p  # Mass of current aggregate
            frict_p = aggs[k].f  # Friction factor of current aggregate
            Kn = aggs[k].Kn  # Knudsen-Number of current aggregate
            beta = aggs[k].beta  # Relaxation time
            t_s = aggs[k].t_s  # Scaling time (calculated from RMSD)
            r_s = aggs[k].r_s  # Scaling length
            rot_attempts = 0  # Rotationsparameter which counts how often the direction got changed
            agg_cross = False  # Crossed the Box Boundaries
            cross = False

            # Reset velocities and displacements
            v_x0, v_y0, v_z0 = 0, 0, 0
            dis_x, dis_y, dis_z = 0, 0, 0

            # Random coordinate in x and y direction within the L*L square
            dis_x = random.uniform(L_x_1, L_x_2)
            dis_y = random.uniform(L_y_1, L_y_2)
            if 2 in params.debugging:
                dis_x = float(params.debug_coord[0])
                dis_y = float(params.debug_coord[1])
            elif 3 in params.debugging:
                def_x_1 = params.define_area / 2
                def_x_2 = -params.define_area / 2
                def_y_1 = params.define_area / 2
                def_y_2 = -params.define_area / 2

                dis_x = random.uniform(def_x_1, def_x_2)
                dis_y = random.uniform(def_y_1, def_y_2)
            # Localized at a certain hight above ground
            dis_z = max_z + Height_above

            i = 0

            ########################################
            ####### Establish neighbor-list ########
            ########################################
            # The k-Dimensional tree of scipy is used:
            # http://docs.scipy.org/doc/scipy/reference/generated/scipy.spatial.cKDTree.html
            # This is a binary tree which sorts the z-Coordinate of the particle. The list
            # is established in one dimension, because the periodic boundary conditions prevent
            # a neighbor recognition across the borders. Note, that the cKDtree is faster than the KDtree, but work
            # in almost the
            neigh_tree = np.zeros((count_PP, 2))
            indx = 0
            if count_PP > 0:
                for agg in range(len(aggs)):
                    if agg < k:
                        for PP in range(aggs[agg].Np):
                            neigh_tree[indx, 0] = aggs[agg].x[PP]
                            neigh_tree[indx, 1] = aggs[agg].y[PP]
                            indx += 1
                tree_before = datetime.datetime.now()
                tree = spatial.cKDTree(neigh_tree)
                tree_after = datetime.datetime.now()
                building_tree_time.append([tree_after - tree_before])
            ########################################
            ####### Particles start moving #########
            ########################################
            # Initialize aggregate inside box
            aggs[k].shift_to(dis_x, dis_y, dis_z)
            aggs[k].checkWall(BOX)
            ###################################
            ##### Check all the particles #####
            ###################################
            particles_hitted = []

            for PP in range(aggs[k].Np):
                # # at first transfer the particles to pt as a working environment
                pt[count_PP + PP].x = aggs[k].x[PP]
                pt[count_PP + PP].y = aggs[k].y[PP]
                pt[count_PP + PP].z = aggs[k].z[PP]
                # count the neighbors which are potentially hit
                # Therefore the previous created tree is searched for any particle within radius
                # 'search_radius' around the particle location
                # If the aggregate is the first one, there is no tree to search within. Therefore an
                # empty list will be created
                if count_PP > 0:
                    # Search radius extension factor
                    sref = 1.01
                    search_radius = (pt[count_PP + PP].r + maximal_radius) * sref
                    neighbors = tree.query_ball_point([pt[count_PP + PP].x, pt[count_PP + PP].y], search_radius)
                    if pt[count_PP + PP].x < (BOX[0] + maximal_radius + pt[count_PP + PP].r):
                        addneighbors = tree.query_ball_point([pt[count_PP + PP].x - 2 * BOX[0], pt[count_PP + PP].y],
                                                             search_radius)
                        for neigh in addneighbors:
                            neighbors.append(neigh)
                    elif pt[count_PP + PP].x > (BOX[1] - maximal_radius - pt[count_PP + PP].r):
                        addneighbors = tree.query_ball_point([pt[count_PP + PP].x - 2 * BOX[1], pt[count_PP + PP].y],
                                                             search_radius)
                        for neigh in addneighbors:
                            neighbors.append(neigh)
                    if pt[count_PP + PP].y < (BOX[2] + maximal_radius + pt[count_PP + PP].r):
                        addneighbors = tree.query_ball_point([pt[count_PP + PP].x, pt[count_PP + PP].y - 2 * BOX[2]],
                                                             search_radius)
                        for neigh in addneighbors:
                            neighbors.append(neigh)
                    elif pt[count_PP + PP].y > (BOX[3] - maximal_radius - pt[count_PP + PP].r):
                        addneighbors = tree.query_ball_point([pt[count_PP + PP].x, pt[count_PP + PP].y - 2 * BOX[3]],
                                                             search_radius)
                        for neigh in addneighbors:
                            neighbors.append(neigh)
                else:
                    neighbors = []
                # the coordinates pt[np_box].x/y/z will be changed by set_new_pos in case of hitting surface or particle
                hit = f.set_new_pos(pt, (count_PP + PP), neighbors, BOX)
                aggs[k].hit[PP] = aggs[k].z[PP] - pt[count_PP + PP].z
            ############### done ###############
            ''' Check if the entire aggregate is outside the box and make the transition
            to the other side '''

            # Check which particle in the agg did the hit (closest to structure in
            # direction of displacement), if no hit the first particle will be remembered
            dummy = 1e99
            print()
            for PP in range(aggs[k].Np):
                if aggs[k].hit[PP] < dummy:
                    dummy = aggs[k].hit[PP]
                    agg_p_hit = PP

            # Give all particles in the agg the position according to the displacement or hit
            for PP in range(aggs[k].Np):
                # If this if loop is not integrated single particles will be moved even though they cross the border.
                # Like this they are not moved back to their position which was changed within the wall_check function.
                aggs[k].z[PP] = aggs[k].z[PP] - aggs[k].hit[agg_p_hit]
                pt[count_PP + PP].z = aggs[k].z[PP]
                pt[count_PP + PP].agg = k
            count_PP += aggs[k].Np
            ''' If hit occured assign all particles of the agg to the structure,
            rank the structure and then find the highest particle '''
            max_z = f.get_max_z(pt, count_PP, max_z)

        # Write the trajectory if it is demanded
        f.write_trajectory(aggs, BOX, params.rmin, max_z, count_steps, k, write_lammpstrj, params.modules,
                           params.Material_selection)
        # write_lammpstrj.close()

    if 9 not in params.debugging and Redirection == 0:
        # Progress bar
        sys.stdout.write(f.progress(k + 1, len(aggs)) + '\r')
        sys.stdout.flush()

    # write_lammpstrj.close()
    return max_z

def post_proc(params, pt, aggs, source_file, max_z):
    output_folder = params.output_folder
    N_PP, N_tot, N_PP_per_Agg = statistics.calc_N(params, source_file, params.N_tot)
    Pe = 0
    Vel = 0
    ###--- Output of the Primary Particle Data (coordinates, radii, agg, etc.)
    # Box params
    L_x_1 = - params.Box_size / 2
    L_x_2 = params.Box_size / 2
    L_y_1 = - params.Box_size / 2
    L_y_2 = params.Box_size / 2
    Height_above = params.Height_above
    BOX = [L_x_1, L_x_2, L_y_1, L_y_2]

    # Standard outputfiles
    now = datetime.datetime.now().replace(microsecond=0)
    current_date = now.strftime("%Y-%m-%d %H:%M")
    rmin, rmax, sigma, r_mean = statistics.calc_r(params, source_file)
    r_ref = r_mean * pow(10, 9)
    r_ref_m = r_ref * pow(10, -9)
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    # Output of the Aggregate and Simulation Data

    output_1_PP = output_folder + params.experiment + '_Particles.dat'
    output_1_AGG = output_folder + params.experiment + '_SimAgg.dat'
    write_1_1 = open(output_1_PP, 'w')
    write_1_1.write('# Here you find just the coordinates of the particles.\n')
    write_1_1.write('# Look into %s for general simulation data an aggregate data\n\n' % (params.experiment + '_SimAgg.dat'))
    write_1_1.write('# Date:\t %s\n' % current_date)
    write_1_1.write('# ID\tAGG\tType\tID_i_Agg\tr (m),x (m),y (m),z (m)\n')
    INDEX = 0
    for PP in range(len(pt)):
        write_1_1.write('%i\t%i\t%i\t%i\t%1.18e\t%1.18e\t%1.18e\t%1.18e\n' \
                        % ((pt[PP].idx + 1), pt[PP].agg, pt[PP].TYPE, pt[PP].pp, \
                           pt[PP].r * r_ref_m, pt[PP].x * r_ref_m, pt[PP].y * r_ref_m, pt[PP].z * r_ref_m))
    write_1_1.close()

    write_1_2 = open(output_1_AGG, 'w')
    write_1_2.write('# Here you find the simulation outputdata\n')
    write_1_2.write('# Look into %s for general simulation data and aggregate data\n#\n' % (params.experiment + '_Particles.dat'))
    write_1_2.write('#######################################################\n')
    write_1_2.write('# Experiment:\t\t\t\t\t\t\t%s\n' % params.experiment)
    write_1_2.write('# Date:\t\t\t\t\t\t\t\t\t%s\n' % current_date)
    write_1_2.write('\n#######################################################\n')
    write_1_2.write('# Conditions acc to:\t\t\t\t\t%s\n' % params.Mov_Method)
    write_1_2.write('# Peclet Number:\t\t\t\t\t\t%1.4e\n' % Pe)
    write_1_2.write('# z-Velocity (m/s):\t\t\t\t\t\t%1.4e\n' % Vel)
    write_1_2.write('# Box_size (rel. to r_mean):\t\t\t%i\n' % params.Box_size)
    write_1_2.write('# Hight above max(z) (rel. to r_mean):\t%i\n' % params.Height_above)
    write_1_2.write('# Timestep adjustment:\t\t\t\t\t%1.2f\n' % params.dt)
    write_1_2.write('# Temperature (K):\t\t\t\t\t\t%3.2f\n' % params.T)
    write_1_2.write('#######################################################\n')
    write_1_2.write('# Particle Size Distribution:\t\t\t%s\n' % params.r_dist)
    if params.r_dist == 'random':
        write_1_2.write('# Particle Minimum Radius (nm):\t\t\t%2.2f\n' % params.rmin)
        write_1_2.write('# Particle Maxyimum Radius (nm):\t\t\t%2.2f\n' % params.rmax)
    elif params.r_dist == 'none':
        write_1_2.write('# Particle Radius (nm):\t\t\t\t\t%2.2f\n#\n' % params.rmin)
    else:
        write_1_2.write('# Particle Mean Radius (nm):\t\t\t%2.2f\n' % params.rmin)
        write_1_2.write('# Particle sDerivation:\t\t\t\t\t%2.2f\n' % params.r_sigma)
    if params.Fric_method == 1:
        write_1_2.write('# Calculation Friction acc to:\t\t\tChen & Dahneke\n')
    elif params.Fric_method == 2:
        write_1_2.write('# Calculation Friction acc to:\t\t\tCunnigham\n')
    elif params.Fric_method == 3:
        write_1_2.write('# Calculation Friction acc to:\t\t\tZhang et al.\n')
    write_1_2.write('#######################################################\n')
    write_1_2.write('# Primary particles:\t\t\t\t\t%i\n' % params.N_tot)
    write_1_2.write('# Aggregates:\t\t\t\t\t\t\t%i\n' % len(aggs))
    write_1_2.write('# Particle Number Distribution:\t\t\t%s\n' % params.N_dist)
    if params.N_dist == 'none':
        write_1_2.write('# Particles per Aggregate:\t\t\t\t%i\n#\n' % params.N_min)
    elif params.N_dist == 'random':
        write_1_2.write('# Minimum Particles per Aggregate:\t\t%i\n' % params.N_min)
        write_1_2.write('# Maximum Particles per Aggregate:\t\t%i\n' % params.N_max)
    else:
        write_1_2.write('# Mean Paricles per Aggregate:\t\t\t%i\n' % params.N_min)
        write_1_2.write('# Standard Derivation:\t\t\t\t\t%i\n' % params.N_sigma)
    write_1_2.write('# Fractal Dimension / Prefactor:\t\t%1.2f\t%1.2f\n' % (params.D_f, params.k_f))
    write_1_2.write(
        '# Epsilon,Delta,Deformation:\t\t\t%1.3f\t%1.3f\t%1.3f\n' % (params.epsilon, params.delta, params.deformation))
    if 4 in params.modules:
        write_1_2.write('# Material selection method:\t\t\t%i\n' % params.Material_selection)
        write_1_2.write('# Density Material 1 (kg/m^3):\t\t\t\t%1.2e\n' % params.rho_1)
        write_1_2.write('# Density Material 2 (kg/m^3):\t\t\t\t%1.2e\n' % params.rho_2)
    else:
        write_1_2.write('# One Material is used\n')
        write_1_2.write('# Density (kg/m^3):\t\t\t\t\t\t%1.2e\n#\n' % params.rho_1)
    write_1_2.write('#######################################################\n')
    write_1_2.write('# Agg\tNp\tDiff_p(m^2/s)\tfrict\tm_p(kg)\tt_s(s)\tr_g\trho(kg/m^3)\n')
    for agg in range(len(aggs)):
        write_1_2.write('%i\t%i\t%1.6e\t%1.6e\t%1.6e\t%1.6e\t%1.6e\t%1.3e\n' % (
        agg, aggs[agg].Np, aggs[agg].D, aggs[agg].f, aggs[agg].m_p, aggs[agg].t_s, aggs[agg].r_g, aggs[agg].rho))
    write_1_2.close()

    ###--- Output of the last state as trajectory file to visualize in ovito
    output_2_Traj = output_folder + params.experiment + '_final' + params.traj_format

    with open(output_2_Traj, 'w') as write_final:
        write_final.write('ITEM: NUMBER OF ATOMS\n%i\n' % params.N_tot)
        write_final.write('ITEM: BOX BOUNDS pp pp pp\n')
        write_final.write('%1.3e %1.3e\n' % (BOX[0] * r_ref, BOX[1] * r_ref))
        write_final.write('%1.3e %1.3e\n' % (BOX[2] * r_ref, BOX[3] * r_ref))
        write_final.write('0 %1.3e\n' % (max_z * r_ref))
        write_final.write('ITEM: ATOMS id type x y z Radius c_tot c_neigh c_same\n')
        INDEX = 0
        for agg in range(len(aggs)):
            for PP in range(aggs[agg].Np):
                write_final.write('%i %i %1.6e %1.6e %1.6e %1.6e %i %i %i\n' % (
                INDEX, agg, aggs[agg].x[PP] * r_ref, aggs[agg].y[PP] * r_ref, aggs[agg].z[PP] * r_ref,
                aggs[agg].r[PP] * r_ref, pt[INDEX].c_0[0], pt[INDEX].c_0[1], pt[INDEX].c_0[2]))
                INDEX += 1
    write_final.close()


    ###--- Paraview VTK file
    if 5 in params.output:
        outputfile = output_folder + params.experiment + '.vtk'
        with open(outputfile,'w') as write:
            write.write('# vtk DataFile Version 3.0\n')
            write.write('Generated by DLA Aggregate Deposition\n')
            write.write('ASCII\n')
            write.write('DATASET POLYDATA\n')
            write.write('POINTS %i float\n' %N_PP)
            for agg in aggs:
                for PP in range(agg.Np):
                    write.write('%1.6e %1.6e %1.6e\n' %(agg.x[PP]*r_ref*pow(10,-9), agg.y[PP]*r_ref*pow(10,-9), agg.z[PP]*r_ref*pow(10,-9)))
            write.write('VERTICES %i %i\n' %(N_PP, 2*N_PP))
            for PP in range(N_PP):
                write.write('1 %i\n' %PP)
            write.write('POINT_DATA %i\n' %N_PP)
            write.write('SCALARS radius float 1\n')
            write.write('LOOKUP_TABLE default\n')
            for agg in aggs:
                for PP in range(agg.Np):
                    write.write('%1.6e\n' %(agg.r[PP]*r_ref*pow(10,-9)))
            write.write('SCALARS type float 1\n')
            write.write('LOOKUP_TABLE default\n')
            for agg in aggs:
                for PP in range(agg.Np):
                    write.write('%i\n' %(agg.TYPE[PP]))
            write.write('SCALARS id float 1\n')
            write.write('LOOKUP_TABLE default\n')
            for PP in range(N_PP):
                write.write('%i\n' %(PP))
            write.write('SCALARS aggregate float 1\n')
            write.write('LOOKUP_TABLE default\n')
            for agg in range(len(aggs)):
                for PP in range(aggs[agg].Np):
                    write.write('%1.1f\n' %(agg))
            write.write('SCALARS c_tot float 1\n')
            write.write('LOOKUP_TABLE default\n')
            for agg in aggs:
                for PP in range(agg.Np):
                    write.write('%i\n' %(agg.c_0[PP,0]))
            write.write('SCALARS c_neig float 1\n')
            write.write('LOOKUP_TABLE default\n')
            for agg in aggs:
                for PP in range(agg.Np):
                    write.write('%i\n' %(agg.c_0[PP,1]))
            write.write('SCALARS c_same float 1\n')
            write.write('LOOKUP_TABLE default\n')
            for agg in aggs:
                for PP in range(agg.Np):
                    write.write('%i\n' %(agg.c_0[PP,2]))
            #if 3 in params.modules:
            #    for path in range(len(percolation_paths)):
            #        write.write('SCALARS perc_%i float 1\n' %(path+1))
            #        write.write('LOOKUP_TABLE default\n')
            #        for PP in range(len(pt)):
            #            if pt[PP].idx in percolation_paths[path]:
            #                write.write('1\n')
            #            else:
            #                write.write('0\n')

            write.close()
        outputfile = output_folder + params.experiment + '_box.vtk'
        with open(outputfile,'w') as write:
            write.write('# vtk DataFile Version 3.0\n')
            write.write('Generated by DLA Aggregate Deposition\n')
            write.write('ASCII\n')
            write.write('DATASET RECTILINEAR_GRID\n')
            write.write('DIMENSIONS 2 2 2\n')
            write.write('X_COORDINATES 2 float\n')
            write.write('%1.2e %1.2e\n' %((BOX[0]*r_ref)*pow(10,-9), (BOX[1]*r_ref)*pow(10,-9)))
            write.write('Y_COORDINATES 2 float\n')
            write.write('%1.2e %1.2e\n' %((BOX[2]*r_ref)*pow(10,-9), (BOX[3]*r_ref)*pow(10,-9)))
            write.write('Z_COORDINATES 2 float\n')
            write.write('0 %1.2e\n' %(max_z*r_ref*pow(10,-9)))
            write.close()

        ###### - Output Henrike
        if 6 in params.output:
            output = output_folder + params.experiment + '_mixing-state' + params.traj_format
            with open(output,'w') as write_henrike:
                write_henrike.write('ITEM: TIMESTEP\n%i\n' %(count_steps))
                write_henrike.write('ITEM: NUMBER OF ATOMS\n%i\n' %params.N_tot)
                write_henrike.write('ITEM: BOX BOUNDS pp pp pp\n')
                write_henrike.write('%1.3e %1.3e\n' %(BOX[0]*r_ref, BOX[1]*r_ref))
                write_henrike.write('%1.3e %1.3e\n' %(BOX[2]*r_ref, BOX[3]*r_ref))
                write_henrike.write('0 %1.3e\n' %(max_z*r_ref))
                write_henrike.write('ITEM: ATOMS id Aggregate Cluster type x y z Radius\n')
                INDEX = 0
                for agg in range(len(aggs)):
                        for PP in range(aggs[agg].Np):
                            write_henrike.write('%i %i %i %i %1.6e %1.6e %1.6e %1.6e\n' %(INDEX,agg,
                                                aggs[agg].cluster[PP],
                                                aggs[agg].TYPE[PP],
                                                aggs[agg].x[PP]*r_ref,
                                                aggs[agg].y[PP]*r_ref,
                                                aggs[agg].z[PP]*r_ref,
                                                aggs[agg].r[PP]*r_ref))
                            INDEX += 1
            write_henrike.close()

        ###### - k-d Tree Time
        if 9 in params.output:
            output = output_folder + params.experiment + '_tree-time.dat'
            with open(output,'w') as write:
                write.write('# Time\n')
                for lin in range(len(building_tree_time)):
                    write.write('%i\t%s\n'%((lin+1),building_tree_time[lin][0]))

        '''
        #####################################################
        ### Remove the copy of the inputfile and close files
        #####################################################
        '''
        if os.getcwd() != source_file:
            if target != source_file + '/' + 'DLA_input.py':
                os.remove(target)




