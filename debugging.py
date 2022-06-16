import numpy as np
import functions as f
import sys, fcntl, termios, struct
from scipy import spatial

if sys.stdout.isatty():
    COLS = struct.unpack('hh',  fcntl.ioctl(sys.stdout, termios.TIOCGWINSZ, '1234'))[1] - 1
else:
    COLS = 60
"""
    Here are some debugging routines
"""
def write_standard_particle(a):
    writing = open('standard_particle.py','w')
    writing.write('N = %i\nD = %1.6e\nm_p = %1.6e\nfrict_p = %1.6e\nKn =%1.6e\nrho = %1.3e\n' %(a.Np,a.D,a.m_p,a.f,a.Kn,a.rho))
    writing.write('r_g = %1.6e\nm = %1.6e\nr_small = %1.6e\nloops = %i\nmax_i = %i\n' %(a.r_g, a.nd, a.r_small, a.loops,a.max_i))
    writing.write('beta = %1.6e\nt_s = %1.6e\nr_s = %1.6e\n' \
        %(a.beta, a.t_s, a.r_s))
    writing.write('r = [%1.6e' %a.r[0])
    for ii in range(len(a.r)-1):
        writing.write(', %1.6e' %a.r[ii+1])
    writing.write(']\n')

    writing.write('x = [%1.6e' %a.x[0])
    for ii in range(len(a.x)-1):
        writing.write(', %1.6e' %a.x[ii+1])
    writing.write(']\n')

    writing.write('y = [%1.6e' %a.y[0])
    for ii in range(len(a.y)-1):
        writing.write(', %1.6e' %a.y[ii+1])
    writing.write(']\n')

    writing.write('z = [%1.6e' %a.z[0])
    for ii in range(len(a.z)-1):
        writing.write(', %1.6e' %a.z[ii+1])
    writing.write(']')
    writing.close()

def check_spaces(pt,BOX):
    errors = []
    epsilon = 0.01
    for PP in range(len(pt)):
        for check_pt in range(len(pt)):
            if (check_pt != PP) and [pt[PP].idx, pt[check_pt].idx] not in errors:
                check_distance = f.distance_pt(pt[check_pt], pt[PP])
                if check_distance < (pt[check_pt].r + pt[PP].r-epsilon):
                      if PP < check_pt:
                        errors.append([int(pt[PP].idx), int(pt[check_pt].idx), float(check_distance)])
    return errors

def check_overlapping(name,output_folder,inputfile,traj_format):
    """
        This scripts imports the _final.lammpstrj file from DLA Simulations and searches for overlapping Particles
    """

    print('Starting Debugging routine: Search for overlappings ...', end='')
    sys.stdout.flush()
    outputfile = output_folder + name + '_overlappings.dat'
    lammps = output_folder + name + '_overlappings'+traj_format

    # --- Functions
    def distance_sq(PP,dummy):
        distance = ( pow(DATA[PP][2]-DATA[dummy][2],2) + pow(DATA[PP][3]-DATA[dummy][3],2) + pow(DATA[PP][4]-DATA[dummy][4],2) )
        return distance

    # --- Open Datafile
    read = open(inputfile)
    lines = read.readlines()
    read.close()
    DATA = []

    for lin in range(9, len(lines)):
        tmp = lines[lin] .split()
        DATA.append([int(tmp[0]), int(tmp[1]), float(tmp[2]), float(tmp[3]), float(tmp[4]), float(tmp[5])])

    ERROR = []
    INDICES = []
    for PP in range(len(DATA)):
        for dummy in range(len(DATA)):
            if PP != dummy:
                dist = round(distance_sq(PP,dummy),1)
                if dist < round(pow(DATA[PP][5] + DATA[dummy][5],2),1)  and \
                [dummy, PP, np.sqrt(dist)] not in ERROR:
                    INDICES.append(PP)
                    INDICES.append(dummy)
                    ERROR.append([PP, dummy, np.sqrt(dist)])

    if len(ERROR) > 0:
        write = open(outputfile,'w')
        write.write('Overlappings found in %s\n' %inputfile)
        write.write('#\tPP1\tPP2\tDistance (m)\n')
        for lin in range(len(ERROR)):
            write.write('%i\t%i\t%i\t%e\n' %(lin+1, ERROR[lin][0],ERROR[lin][1],ERROR[lin][2]))
        write.close()

        write = open(lammps,'w')
        write.write('ITEM: TIMESTEP\n0\n')
        write.write('ITEM: NUMBER OF ATOMS:\n%i\n' %len(DATA))
        write.write('ITEM: BOX BOUNDS pp pp pp\n')
        write.write(lines[5])
        write.write(lines[6])
        write.write(lines[7])
        write.write('ITEM: ATOMS id type x y z Radius\n')
        for lin in range(len(DATA)):
            if lin in INDICES:
                write.write('%i 1 %e %e %e %e\n' %(lin, DATA[lin][2],DATA[lin][3],DATA[lin][4],DATA[lin][5]))
            else:
                write.write('%i 0 %e %e %e 1\n' %(lin, DATA[lin][2],DATA[lin][3],DATA[lin][4],DATA[lin][5]))
        write.close()
    print(' done\n\nFound %i errors' %(len(ERROR)))
    return len(ERROR)


def check_boundaries(name,output_folder,inputfile, N_tot,traj_format):
    ''' Multiply the system to check boundary overlappings'''

    print('\n%s' %f.delim_char(COLS,'*'))
    print('Search for overlappings including boundaries\n')
    outputfile = output_folder + name + '_boundary_overlappings.dat'
    outputfile_name = name + '_boundary_overlappings.dat'
    lammps = output_folder + name + '_boundary_overlappings'+traj_format
    lammps_name = name + '_boundary_overlappings'+traj_format

    # --- Open Datafile
    read = open(inputfile)
    lines = read.readlines()
    read.close()
    data = []
    tree_data = []

    def distance(PP,dummy):
        distance = np.sqrt( pow( data[PP][1]-data[dummy][1] ,2) + pow( data[PP][2]-data[dummy][2] ,2) + pow( data[PP][3]-data[dummy][3] ,2) )
        return distance

    def allowed_distance(PP,dummy):
        allowed_distance = data[PP][4]+data[dummy][4]
        return allowed_distance

    BOX_X = float(lines[5].split()[-1]) * 2
    BOX_Y = float(lines[6].split()[-1]) * 2
    max_r = 0
    for lin in range(9, len(lines)):
        tmp = lines[lin].split()
        data.append([int(tmp[1]), float(tmp[2]), float(tmp[3]), float(tmp[4]), float(tmp[5])])
        tree_data.append([float(tmp[1]),float(tmp[2]),float(tmp[3])])
        if data[-1][4] > max_r:
            max_r = data[-1][4]
    particles = len(data)
    for dat in range(len(data)):
        data.append([data[dat][0],data[dat][1]+BOX_X,data[dat][2],data[dat][3],data[dat][4]])
        data.append([data[dat][0],data[dat][1],data[dat][2]+BOX_Y,data[dat][3],data[dat][4]])
        tree_data.append([data[dat][1]+BOX_X,data[dat][2],data[dat][3]])
        tree_data.append([data[dat][1],data[dat][2]+BOX_Y,data[dat][3]])

    tree = spatial.cKDTree(tree_data)

    total_overlappings = 0
    overlapping_list = []
    ERROR_list = []
    INDICES = []
    count = 1
    for PP in range(particles):
        neighbor = tree.query_ball_point(tree_data[PP],data[PP][4]+max_r)
        overlappings = 0
        for part in neighbor:
            if part != PP and [part, PP] not in overlapping_list:
                if distance(PP,part) < 0.99 * allowed_distance(PP,part):
                    total_overlappings += 1
                    overlappings += 1
                    overlapping_list.append([PP, part])
                    ERROR_list.append([count, PP, part,distance(PP,part),data[PP][4],data[part][4]])
                    INDICES.append(PP)
                    INDICES.append(part)
                    count += 1


    write = open(outputfile,'w')
    write.write('Overlappings found in %s\n' %inputfile)
    write.write('#\tPP1\tPP2\tDistance (m)\tRadius PP1\tRadius PP2\n')
    for lin in range(len(ERROR_list)):
        write.write('%i\t%i\t%i\t%e\t%e\t%e\n' %(lin+1, ERROR_list[lin][1],ERROR_list[lin][2],ERROR_list[lin][3],ERROR_list[lin][4],ERROR_list[lin][5], ))
    write.close()

    write = open(lammps,'w')
    write.write('ITEM: TIMESTEP\n0\n')
    write.write('ITEM: NUMBER OF ATOMS:\n%i\n' %len(data))
    write.write('ITEM: BOX BOUNDS pp pp pp\n')
    write.write(lines[5])
    write.write(lines[6])
    write.write(lines[7])
    write.write('ITEM: ATOMS id type x y z Radius\n')
    for lin in range(len(data)):
        if lin in INDICES:
            write.write('%i 1 %e %e %e %e\n' %(lin, data[lin][1],data[lin][2],data[lin][3],data[lin][4]))
        else:
            write.write('%i 0 %e %e %e 1\n' %(lin, data[lin][1],data[lin][2],data[lin][3]))
    write.close()
    print('>>> %s errors found' %len(ERROR_list))
    print('\n%s\n' %f.delim_char(COLS,'*'))
    return len(ERROR_list)
