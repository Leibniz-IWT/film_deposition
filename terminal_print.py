import functions as f
import sys, fcntl, termios, struct

if sys.stdout.isatty():
    Redirection = 0
    COLS = struct.unpack('hh',  fcntl.ioctl(sys.stdout, termios.TIOCGWINSZ, '1234'))[1] -1
else:
    COLS = 60
    Redirection = 1


def loghead(params, name, PID, num_cores, now, N_PP, N_PP_per_Agg,rmin,rmax):
    ##### General Output for the logs #####
    print('\n%s'%f.delim_char(COLS,'*'))
    print('Process data')
    print('\tExperiment:\t\t%s' %name)
    print('\tStarting Time:\t\t%s' %now.strftime("%Y-%m-%d %H:%M"))
    print('\tPID:\t\t\t%s'%PID)
    print('\tCPUs:\t\t\t%s'%num_cores)
    print('\n%s\nSystem parameters' %f.delim_char(COLS,'*'))
    if params.Mov_Method == 'Pe':
        print('\tPe:\t\t\t%1.2e' %params.Pe)
    elif params.Mov_Method == 'Vel':
        print('\tv_z:\t\t\t-%1.2e m/s' %params.Vel)
    print('\tBox size:\t\t%i' %params.Box_size)
    print('\tPrimary Particles:\t%i' %N_PP)
    print('\tCluster per Aggregate:\t%i' %params.Cluster)
    print('\tAggregates:\t\t%i' %len(N_PP_per_Agg))
    print('\tParticle Number Dist.:\t%s' %params.N_dist)
    if params.N_dist == 'predef':
        print('\t\tfile:\t\t%s' %params.pre_N_file)
    print('\tParticle Size Dist.:\t%s' %params.r_dist)
    if params.r_dist == 'norm':
        print('\tMedian radius (nm):\t%1.3e' %rmin)
        print('\tSigma (nm):\t\t%1.3e' %params.r_sigma)
    elif params.r_dist == 'lognorm':
        print('\tMedian radius (nm):\t%1.3e' %rmin)
        print('\tSigma:\t\t\t%1.3e' %params.r_sigma)
    elif params.r_dist == 'random':
        print('\trmin (m):\t\t%1.3e' %rmin)
        print('\trmax (m):\t\t%1.3e' %rmax)
    elif params.r_dist == 'predef':
        print('\trmin (m):\t\t%1.3e' %rmin)
        print('\t\tfile:\t\t%s' %params.pre_dp_file.split('/')[-1])
    else:
        print('\trmin (m):\t\t%1.3e' %rmin)
        print('\trmax (m):\t\tnone')
    print('\tNeighbor Mode:\t\t%i' %params.coord_threshold_mode)
    if params.coord_threshold_mode == 3:
        print('\tThreshold (nm):\t\t%i' %params.coord_threshold)
    if params.Fric_method == 1:
        print('\tFriction Model:\tChan&Dahneke')
    elif params.Fric_method == 2:
        print('\tFriction Model:\tCunningham')
    elif params.Fric_method == 3:
        print('\tFriction Model:\t\tZhang (Tabulated)')
        print('\t\tfile\t\t%s' %params.diffusionfile)




    print('\n%s' %f.delim_char(COLS,'*'))
    print('Modules')
    if not params.modules or (len(params.modules) == 1 and 0 in params.modules):
        print('\tNo modules selected')
    if 1 in params.modules:
        print('\tMultiple Deposition (1)')
    if 2 in params.modules:
        print('\tCalculate Coordination (2)')
        print('\t\tThreshold: ', end='')
        if params.coord_threshold_mode == 1:
            print('One rp (1)')
        elif params.coord_threshold_mode == 2:
            print('Delta (2)')
        elif params.coord_threshold_mode == 3:
            print('Custom: %f nm (3)' %params.coord_threshold)
    if 3 in params.modules:
        print('\tCalculate Percolation (3)')
        print('\t\tThreshold: ', end='')
        if params.perc_threshold_mode == 1:
            print('One rp (1)')
        elif params.perc_threshold_mode == 2:
            print('Delta (2)')
        elif params.perc_threshold_mode == 3:
            print('Custom: %f nm (3)' %params.perc_threshold)
    if 4 in params.modules:
        print('\tMulticomponent system (4)')
    print('\n%s' %f.delim_char(COLS,'*'))
    print('Outputfiles')
    print('\tDefault outputdata')
    if 1 in params.output:
        print('\tlog file (1)')
    if 2 in params.output:
        print('\tLIGGGHTS hybrid granular/molecular data file (2)')
    if 3 in params.output:
        print('\tLIGGGHTS granular data file (3)')
    if 4 in params.output:
        print('\tWhole trajectory (stepsize %i) (4)' %params.steps_trajectory)
    if 5 in params.output:
        print('\tVTK Output for paraview (5)')
    if 6 in params.output:
        print('\tMixing State File (6)')
    if 7 in params.output:
        print('\tHistogram with primary particle size distibution (7)')
    if 8 in params.output:
        print('\tHistogram with primary particle number distibution (8)')
    if 9 in params.output:
        print('\tBuilding Time of kd-tree (9)')
    if 10 in params.output:
        print('\tLAMMPS data file molecular (10)')
    print('\n%s' %f.delim_char(COLS,'*'))
    print('Debugging routines:')
    if not params.debugging or (len(params.debugging) == 1 and 0 in params.debugging):
        print('\tnone')
    else:
        if 1 in params.debugging:
            print('\tOnly aggregate generation (1)')
        if 2 in params.debugging:
            print('\tStart at coordinate [%3.1f, %3.1f] (2)' %(params.debug_coord[0],params.debug_coord[1]))
        if 3 in params.debugging:
            print('\tStart in a defined area x,y = +- %3.1f (3)' %(params.define_area/2))
        if 4 in params.debugging:
            print('\tCalculate overlapping particles with boundaries (4)')
        if 5 in params.debugging:
            print('\tCheck overlapping after each deposition (5)')
        if 6 in params.debugging:
            print('\tQuit in case of overlapping (6)')
        if 7 in params.debugging:
            print('\tUse a custom seed for the trajectory (%i) (7)' %params.custom_seed_traj)
        if 8 in params.debugging:
            print('\tUse a custom seed for the aggregates (%i) (8)' %params.custom_seed_aggs)
        if 9 in params.debugging:
            print('\tSuppress progress bar (9)')
    print('\n%s' %f.delim_char(COLS,'*'))
    print('\n##### Aggregate generation initiated')
    sys.stdout.flush()
