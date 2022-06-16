import numpy as np
import random,sys,fcntl,termios,struct,time,matplotlib,os,smtplib
import matplotlib.pyplot as plt
from scipy import spatial

""" Here you'll find functions which are used throughout the code"""

#############################################
########## Distance of two points ###########
#############################################

def distance(x2,x1,y2,y1,z2,z1):
    distance = np.sqrt(pow((x2-x1),2)+pow((y2-y1),2)+pow((z2-z1),2))
    return distance

# The squared distance calculation is faster than using the sqareroot!
def distance_sq(x2,x1,y2,y1,z2,z1):
    distance_sq = (pow((x2-x1),2)+pow((y2-y1),2)+pow((z2-z1),2))
    return distance_sq

def distance_pt(pt_2,pt_1):
    x2, x1 = pt_2.x, pt_1.x
    y2, y1 = pt_2.y, pt_1.y
    z2, z1 = pt_2.z, pt_1.z
    distance = np.sqrt(pow((x2-x1),2)+pow((y2-y1),2)+pow((z2-z1),2))
    return distance

def distance_pt_sq(pt_2,pt_1):
    x2, x1 = pt_2.x, pt_1.x
    y2, y1 = pt_2.y, pt_1.y
    z2, z1 = pt_2.z, pt_1.z
    distance = (pow((x2-x1),2)+pow((y2-y1),2)+pow((z2-z1),2))
    return distance
#############################################
######## Random Standard Derivation##########
#############################################


""" A function to calculate a random standard derivation for gauss distributions """
def rnd_gauss_deriv(params):
    method = params.sdev_method
    # Box, G.E.P, M.E. Muller 1958; A note on the generation of random normal deviates, Annals Math. Stat, V. 29, pp. 610-611
    if method == 1:
        r1 = random.random()
        r2 = random.random()
        x1 = np.sqrt(-2 * np.log(r1)) * np.cos(2*np.pi*r2)
        x2 = np.sqrt(-2 * np.log(r1)) * np.sin(2*np.pi*r2)
        return x1
    elif method == 2:
    # A method used by Mädler 2006
        w = 1
        while np.sqrt(w) >= 1:
            x1 = 2 * random.random() - 1
            x2 = 2 * random.random() - 1
            w = x1*x1 + x2*x2
        w = np.sqrt ((-2 * np.log(w))/2)
        y1 = x1 * w
        y2 = x2 * w
        return y1

#############################################
########### Refresh coordinates #############
#############################################

def ref_pt_pos(pt,aggs):
    """ refresh the position of the particles (in pt) after moving an aggregate (in aggs)"""
    for PP in range(len(pt)):
        pt[PP].x = aggs[pt[PP].agg].x[pt[PP].pp]
        pt[PP].y = aggs[pt[PP].agg].y[pt[PP].pp]
        pt[PP].z = aggs[pt[PP].agg].z[pt[PP].pp]

def ref_agg_pos(pt,aggs):
    """ refresh the position of the particles (in agg) after moving (in pt)"""
    PP = 0
    for agg in range(len(aggs)):
        for npp in range(aggs[agg].Np):
            pt[PP].x = aggs[agg].x[npp]
            pt[PP].y = aggs[agg].y[npp]
            pt[PP].z = aggs[agg].z[npp]
            PP += 1

#############################################
################ Get max_z ##################
#############################################

''' Get the highest particle '''
def get_max_z(pt,particles,max_z):
    for PP in range(particles):
        if pt[PP].z > max_z:
            max_z = pt[PP].z
    return max_z

#############################################
################ Check Hit ##################
#############################################

''' Function to search for a hit with another aggregate '''
def set_new_pos(pt, np_set, neighbors,BOX):
    hit = False
    hitDistanceDummy = 1e7
    """
    At first a potential hit with all neighboring particles has to be checked
    This has to be done across the boundaries
    """
    for j in neighbors:
        check_particle = False
        drx = pt[j].r + pt[np_set].r
        dry = pt[j].r + pt[np_set].r
        if abs(pt[np_set].x - pt[j].x) < drx and \
           abs(pt[np_set].y - pt[j].y) < dry:
           check_particle = True
           checkx = pt[j].x
           checky = pt[j].y
        elif abs(pt[np_set].x - (pt[j].x + 2 * BOX[0])) < drx and \
             abs(pt[np_set].y - pt[j].y) < dry:
           check_particle = True
           checkx = pt[j].x + 2*BOX[0]
           checky = pt[j].y
        elif abs(pt[np_set].x - (pt[j].x + 2 * BOX[1])) < drx and \
             abs(pt[np_set].y - pt[j].y) < dry:
           check_particle = True
           checkx = pt[j].x + 2*BOX[1]
           checky = pt[j].y
        elif abs(pt[np_set].x - pt[j].x) < drx and \
             abs(pt[np_set].y - (pt[j].y + 2 * BOX[2])) < dry:
           check_particle = True
           checkx = pt[j].x
           checky = pt[j].y + 2*BOX[2]
        elif abs(pt[np_set].x - pt[j].x) < drx and \
             abs(pt[np_set].y - (pt[j].y + 2 * BOX[3])) < dry:
           check_particle = True
           checkx = pt[j].x
           checky = pt[j].y + 2*BOX[3]
        elif abs(pt[np_set].x - (pt[j].x + 2 * BOX[0])) < drx and \
             abs(pt[np_set].y - (pt[j].y + 2 * BOX[2])) < dry:
           check_particle = True
           checkx = pt[j].x + 2*BOX[0]
           checky = pt[j].y + 2*BOX[2]
        elif abs(pt[np_set].x - (pt[j].x + 2 *BOX[1])) < drx and \
             abs(pt[np_set].y - (pt[j].y + 2 *BOX[2])) < dry:
           check_particle = True
           checkx = pt[j].x + 2*BOX[1]
           checky = pt[j].y + 2*BOX[2]
        elif abs(pt[np_set].x - (pt[j].x + 2 * BOX[0])) < drx and \
             abs(pt[np_set].y - (pt[j].y + 2 * BOX[3])) < dry:
           check_particle = True
           checkx = pt[j].x + 2*BOX[0]
           checky = pt[j].y + 2*BOX[3]
        elif abs(pt[np_set].x - (pt[j].x + 2 * BOX[1])) < drx and \
             abs(pt[np_set].y - (pt[j].y + 2 * BOX[3])) < dry:
           check_particle = True
           checkx = pt[j].x + 2*BOX[1]
           checky = pt[j].y + 2*BOX[3]

        if check_particle:
            ddh = np.sqrt(pow(checkx - pt[np_set].x, 2) + pow(checky - pt[np_set].y, 2))
            hitcondition = ddh <= (pt[j].r + pt[np_set].r)
            if(hitcondition):
                # Remember there is a hit
                hit = True

                ddr = pt[j].r + pt[np_set].r
                dz = pt[j].z + np.sqrt(ddr**2 - ddh**2)
                check_dist = pt[np_set].z - dz

                '''if particle j is the first one, particle np_set hits or
                its traveled distance (trav_dist) towards the hit is shorter
                than to the particle hitted before its new potential position
                is remembered'''
                if(check_dist < hitDistanceDummy):
                    hitDistanceDummy = check_dist

    ''' After checking all particles, the final position of particle np_set is set
    (and coordinates are translated back) next to the closest particle '''
    if(hit == True):
        pt[np_set].z = pt[np_set].z - hitDistanceDummy

    """ Check if the particle hit the surface without hitting a particle"""
    if(hit == False):
        hit = True
        pt[np_set].z = pt[np_set].r
    return hit



def shift_system(x0,y0,z0,BOX,aggs,pt,k):
    for agg in range(len(aggs)):
        if agg < k:
            aggs[agg].store_coords()
    for agg in range(len(aggs)):
        if agg < k:
            for PP in range(aggs[k].Np):
                aggs[agg].x[PP] += x0
                aggs[agg].y[PP] += y0
                aggs[agg].z[PP] += z0

    for agg in range(len(aggs)):
        if agg < k:
            for PP in range(aggs[k].Np):
                if aggs[agg].x[PP] < BOX[0]:
                    aggs[agg].x[PP] -= BOX[0]
                if aggs[agg].x[PP] > BOX[1]:
                    aggs[agg].x[PP] -= BOX[1]
                if aggs[agg].y[PP] < BOX[2]:
                    aggs[agg].y[PP] -= BOX[2]
                if aggs[agg].y[PP] > BOX[3]:
                    aggs[agg].y[PP] -= BOX[3]
    ref_pt_pos(pt,aggs)

def restore_system(aggs,pt,k):
    for agg in range(len(aggs)):
        if agg < k:
            aggs[agg].load_coords()
    ref_pt_pos(pt,aggs)

##########################################################
################## Packing Density #######################
##########################################################

'''
    This function calculates the packing density of 90% of the film height
'''
def packing_density(pt,max_z,lower,upper,BOX):
    LOWER = []
    RELEVANT = []
    if(max_z > upper):
        for PP in range(len(pt)):
            if pt[PP].z < lower:
                LOWER.append([pt[PP].idx,pt[PP].r])
            if pt[PP].z < upper and [pt[PP].idx,pt[PP].r] not in LOWER:
                RELEVANT.append([pt[PP].idx,pt[PP].r])
        VOLUME = 0
        for PP in range(len(RELEVANT)):
            VOLUME += np.pi * 4/3 * pow(pt[PP].r,3)
        rho_pack = VOLUME / ((upper-lower-1)*(2*BOX[1])*(2*BOX[3]))
    else:
        rho_pack = 0
    return rho_pack


##########################################################
##################### Coordination #######################
##########################################################
'''
    This function calculates the coordination number for each particle. Therefore the BOX has to be multiplied to govern periodic boundary conditions
'''
def coordination(pt,aggs,BOX,threshold_mode,threshold,rmin,delta,N_PP,img_name,outputfolder,img_format):
    # box_x = [0,2*BOX[0],2*BOX[1],0,0]
    # box_y = [0,0,0,2*BOX[2],2*BOX[3]]
    tree_data = []
    max_r = 0
    for PP in range(len(pt)):
        tree_data.append([pt[PP].x,pt[PP].y,pt[PP].z])
        if pt[PP].r > max_r: max_r = pt[PP].r
    tree = spatial.cKDTree(tree_data)
    for PP in range(len(pt)):
        neighbors = tree.query_ball_point([pt[PP].x,pt[PP].y,pt[PP].z],pt[PP].r+max_r)
        for dummy in neighbors:
            if dummy != PP:
                # for mirror in range(len(box_x)):
                    # dist = distance(pt[PP].x,pt[dummy].x+box_x[mirror],pt[PP].y,pt[dummy].y+box_y[mirror],pt[PP].z,pt[dummy].z)
                dist = distance_sq(pt[PP].x,pt[dummy].x,pt[PP].y,pt[dummy].y,pt[PP].z,pt[dummy].z)

                '''
                    Three different options to set the threshold:
                    1: dist <=  r1+r2 + r_mean(1,2)
                    2: dist <=  (r1+r2) * delta
                    3  dist <=  r1+r2 + r_threshold
                '''
                if threshold_mode == 1:
                    r_comp = pow(pt[PP].r + pt[dummy].r + (pt[PP].r + pt[dummy].r) / 2,2)
                elif threshold_mode == 2:
                    r_comp = pow((pt[PP].r + pt[dummy].r)*delta,2)
                elif threshold_mode == 3:
                    r_comp = pow((pt[PP].r + pt[dummy].r) + threshold/rmin,2)

                if dist <= r_comp:
                    if pt[PP].agg == pt[dummy].agg:
                        pt[PP].c_0[0] += 1
                        pt[PP].c_0[2] += 1
                    else:
                        pt[PP].c_0[0] += 1
                        pt[PP].c_0[1] += 1

    ''' Calculation of coordination statistics
        coordination[0,0] = total number of contacts in the system
        coordination[0,1] = mean number of contacts per particle
        coordination[1,0] = total number of aggregate-aggregate contacts
        coordination[1,1] = mean number of agg-agg contacts per aggregate
        coordination[2,0] = total number of inner agg contacts
        coordination[2,1] = mean number of contacts within aggregates
    '''
    coordination = np.zeros((3,2))
    x = 0
    for PP in range(len(pt)):
        coordination[0,0] += pt[PP].c_0[0]
        coordination[1,0] += pt[PP].c_0[1]
        coordination[2,0] += pt[PP].c_0[2]
        x += 1

    coordination[0,1] = coordination[0,0]/len(pt)
    coordination[1,1] = coordination[1,0]/len(aggs)
    coordination[2,1] = coordination[2,0]/len(aggs)

    outputfile = outputfolder + img_name + '_coordination.dat'
    with open(outputfile,'w') as write:
        write.write('# File governing the statistics of the coordination\n')
        write.write('Experiment: %s\n' %img_name)
        write.write('Particles: %i\n' %N_PP)
        write.write('Aggregates: %i\n\n' %(len(aggs)))
        write.write('*******************************************************\n')
        write.write('# The calculations include periodic boundary conditions\n')
        write.write('# Therefore the system was mirrored by -x,+x,-y,+y in order to \n')
        write.write('# investigate split aggregates\n')
        write.write('*******************************************************\n\n')
        write.write('Overall particle contacts:\t%i\n' %(coordination[0,0]))
        write.write('Particle contacts per particle:\t%1.2f\n' %(coordination[0,1]))
        write.write('Overall Agg-Agg contacts:\t%i\n' %coordination[1,0])
        write.write('Agg-Agg contacts per agg:\t%1.2f\n' %coordination[1,1])
        write.write('Mean contacts within an agg:\t%1.2f\n' %coordination[2,1])
        write.write('\nPP\tAgg\tc_tot\tc_neigh\tc_sameagg\n')
        for PP in range(len(pt)):
            write.write('%i\t%i\t%i\t%i\t%i\n' %(pt[PP].idx,pt[PP].agg,pt[PP].c_0[0],pt[PP].c_0[2],pt[PP].c_0[1]))
        write.close()


    # Histogram of the coordination
    c_tot = np.zeros(len(pt))
    c_0 = np.zeros(len(pt))
    c_1 = np.zeros(len(pt))
    c_tot_cum = np.zeros(len(pt))
    for PP in range(len(pt)):
        c_tot[PP] = pt[PP].c_0[0]
        c_0[PP] = pt[PP].c_0[2]
        c_1[PP] = pt[PP].c_0[1]
        c_tot_cum[PP] = pt[PP].c_0[1] + pt[PP].c_0[2]
    c = [c_tot,c_1,c_0]
    name = ['Total neighbors', 'Different aggregate', 'Same aggregate']
    color = ['k','r','b']

    f,(ax1,ax2,ax3) = plt.subplots(3,sharex=True)
    co = 0
    for ax in (ax1,ax2,ax3):
        count, bin,ignored = ax.hist(c[co],bins =np.arange(0,6),align='mid',color=color[co], width = 0.7)
        #ax.plot(bin[:-1],count)
        ax.set_title('%s' %name[co])
        ax.set_xlim(0,6)
        ax.set_ylabel('Fraction / -')
        co += 1
    ax3.set_xlabel('Coordination / -')

    count_cum,bin_cum = np.histogram(c_tot_cum,bins =np.arange(0,6))
    ax1.plot(bin_cum[:-1],count_cum,'xr')

    outputfile = outputfolder + img_name + '_coordination' + img_format
    plt.savefig(outputfile)
    return coordination


##########################################################
###################### Percolation #######################
##########################################################

''' Find percolation chains through the film, which connect the top to the bottom '''

def percolation(pt,TWL,delta,threshold_mode,threshold,rmin,name,outputfolder):

    # TWL:              Top_Wall_Limit - The threshold until the chain has to go
    # BWL:              Bottom_Wall_Limit - Should always be zero (bottom of box), unless not the whole system is investigated
    # TWC:              Top_Wall_Contacts - List containing the connecting particles
    # BWC:              Bottom_Wall_Contacts - List containing the connecting particles
    # pt.adj_num:       Number of direct neighbors
    # pt.adj:           List of the indices of the direct neighbors
    # pt.wall_contact:  Is there a wall contact (-1 bottom / +1 top)
    # pt.idx:           Global index
    # pt.agg:           Aggregate index
    # pt.path_chain     Peculation chain
    # pt.predec         Predecessor in the chain

    # Lists and variables
    BWL = 0
    TWC = []
    BWC = []
    path = []
    np_tot = len(pt)
    INF = 1e10
    COMPLETE_PATH = []

    tree_data = []
    max_r = 0
    for PP in range(len(pt)):
        tree_data.append([pt[PP].x,pt[PP].y,pt[PP].z])
        if pt[PP].r > max_r: max_r = pt[PP].r
    tree = spatial.cKDTree(tree_data)

    # Find all direct neighbors for all particles
    for PP in range(len(pt)):
        neighbors = tree.query_ball_point([pt[PP].x,pt[PP].y,pt[PP].z],pt[PP].r+max_r)
        for dummy in neighbors:
            if dummy != PP:
                dist = distance_sq(pt[PP].x,pt[dummy].x,pt[PP].y,pt[dummy].y,pt[PP].z,pt[dummy].z)
                if threshold_mode == 1:
                    r_comp = pow(pt[PP].r + pt[dummy].r + (pt[PP].r + pt[dummy].r) / 2,2)
                elif threshold_mode == 2:
                    r_comp = pow((pt[PP].r + pt[dummy].r)*delta,2)
                elif threshold_mode == 3:
                    r_comp = pow(pt[PP].r + pt[dummy].r + threshold / rmin,2)
                if dist <= r_comp:
                    pt[PP].adj[pt[PP].adj_num] = dummy
                    pt[PP].adj_num += 1
                    pt[dummy].adj[pt[dummy].adj_num] = PP
                    pt[dummy].adj_num += 1

    # Find wall contacts
    for PP in range(np_tot):
        # Top
        if (pt[PP].z >= TWL - pt[PP].r):
            pt[PP].wall_contact = 1
            if not TWC:
                TWC.append(pt[PP].idx)
            if TWC:
                if pt[PP].agg != pt[TWC[-1]].agg:
                    TWC.append(pt[PP].idx)
        # Bottom
        elif (pt[PP].z <= BWL + pt[PP].r):
            pt[PP].wall_contact = -1
            if not BWC:
                BWC.append(pt[PP].idx)
            if BWC:
                if pt[PP].agg != pt[BWC[-1]].agg:
                    BWC.append(pt[PP].idx)


    ''' OUTPUT FILES '''
    write_out = open(outputfolder + name + '_perc.out','w')
    write_out.write('#Agg\tWall_Contact\tPP\n')
    for i in range(np_tot):
        write_out.write('%5i\t%5i\t%5i\t' %(pt[i].agg,pt[i].wall_contact,pt[i].idx))
        for j in range(pt[i].adj_num):
            write_out.write('%5i\t' %pt[i].adj[j])
        write_out.write('\n')
    write_out.close()

    write_path = open(outputfolder + name +'_path.dat','w')
    write_path.write('Path File with format:\n')
    write_path.write('1. Wall_contacts\n2. Index_in_BWC\n3. Global_PP_Index\n###  In case of percolation\n4. Percolation length\n5. Last Particle index\n6. First Particle index\n7. Involving Aggregates\n')
    write_path.write('8. count_path (newline)\n### In case of double break: Path complete\n\n')
    write_path.write('*************************************************************************\n\n')
    write_path.write('##### BOTTOM WALL #####\n')
    ''' OUTPUT FILES '''


    '''
        #################
        Search algorithm
        #################
    '''
    ''' Starting at the Bottom '''
    count_path = 0
    for i in range(len(BWC)):
        write_path.write('%5i\t%5i\t%5i\t' %(len(BWC), i, BWC[0]))
        for k in range(np_tot):
            pt[k].path_chain = 0
            pt[k].dist = INF
            pt[k].predec = INF

        path = []
        path_trj = []
        path.append(BWC[0])
        final_pp = BWC[0]
        pt[BWC[0]].path_chain = 1
        pt[BWC[0]].dist = 0
        BWC.append(BWC[0])
        BWC.pop(0)
        while path and pt[final_pp].wall_contact != 1:
            for j in range(pt[path[0]].adj_num):
                if pt[int(pt[path[0]].adj[j])].path_chain == 0:
                    pt[pt[path[0]].adj[j]].path_chain = 1
                    pt[pt[path[0]].adj[j]].dist = pt[path[0]].dist + 1
                    pt[pt[path[0]].adj[j]].predec = path[0]
                    path.append(pt[pt[path[0]].adj[j]].idx)
                    # End if
                # End for j
            pt[path[0]].path_chain = 2
            final_pp = int(path[0])
            path.pop(0)
            # End while
        # Count aggregates in the chain
        if pt[final_pp].wall_contact == 1:
            count_path += 1
            count_agg = 0
            current_pp = final_pp
            write_path.write('%5i\t' %pt[current_pp].dist)
            write_path.write('%5i\t' %pt[current_pp].idx)

            # Write out the percolation path
            while pt[current_pp].dist > 0:
                current_pp = pt[current_pp].predec
                if pt[current_pp].predec < INF:
                    path_trj.append(pt[current_pp].idx)
                    if pt[current_pp].agg != pt[pt[current_pp].predec].agg:
                        count_agg += 1
                        # End if
                    # End if
                # End while
            COMPLETE_PATH.append(path_trj)
            write_path.write('%5i\t%5i\t\n' %(pt[current_pp].idx, count_agg+1))
            # End if
        # End for i
        write_path.write('\n%5i\n\n' %count_path)

    write_path.write('##### TOP WALL #####\n')
    ''' Top Wall '''
    count_path = 0
    for i in range(len(TWC)):
        write_path.write('%5i\t%5i\t%5i\t' %(len(TWC), i, TWC[0]))
        for k in range(np_tot):
            pt[k].path_chain = 0
            pt[k].dist = INF
            pt[k].predec = INF
            # End for
        path = []
        path_trj = []
        path.append(TWC[0])
        final_pp = int(TWC[0])
        pt[TWC[0]].path_chain = 1
        pt[TWC[0]].dist = 0
        TWC.append(TWC[0])
        TWC.pop(0)
        while path and pt[final_pp].wall_contact != -1:
            for j in range(pt[path[0]].adj_num):
                if pt[pt[path[0]].adj[j]].path_chain == 0:
                    pt[pt[path[0]].adj[j]].path_chain = 1
                    pt[pt[path[0]].adj[j]].dist = pt[path[0]].dist + 1
                    pt[pt[path[0]].adj[j]].predec = path[0]
                    path.append(pt[pt[path[0]].adj[j]].idx)
                    # End if
                # End for j
            pt[path[0]].path_chain = 2
            final_pp = int(path[0])
            path.pop(0)
            # End while
        if pt[final_pp].wall_contact == -1:
            count_path += 1
            count_agg = 0
            current_pp = final_pp
            write_path.write('%5i\t' %pt[current_pp].dist)
            write_path.write('%5i\t' %pt[current_pp].idx)
            while pt[current_pp].dist > 0:
                current_pp = pt[current_pp].predec
                if pt[current_pp].predec < INF:
                    path_trj.append(pt[current_pp].idx)
                    if pt[current_pp].agg != pt[pt[current_pp].predec].agg:
                        count_agg += 1
                        # End if
                    # End if
                # End while
            COMPLETE_PATH.append(path_trj)
            write_path.write('%5i\t%5i\t\n' %(pt[current_pp].idx, count_agg+1))
            # End if
        # End for i
    write_path.write('\n%5i\n\n' %count_path)
    write_path.close()
    return COMPLETE_PATH


##########################################################
##################### Progress Bar #######################
##########################################################
if sys.stdout.isatty():
    COLS = struct.unpack('hh',  fcntl.ioctl(sys.stdout, termios.TIOCGWINSZ, '1234'))[1]

def bold(msg):
    return (u"\033[1m%s\033[0m" %msg)

def progress(current, total):
    prefix = ' Progress: %5i / %5i (%5i %%)' % (current, total, int((current/total)*100))
    bar_start = ' ['
    bar_end = '] '
    if sys.stdout.isatty():
        bar_size = COLS - len(prefix + bar_start + bar_end)
    else:
        bar_size = 10
    amount = int(current / (total / float(bar_size)))
    remain = bar_size - amount
    bar = 'X' * amount + ' ' * remain

    return bold(prefix) + bar_start + bar + bar_end

##########################################################
################### Delimiter strings ####################
##########################################################
def delim_char(length,char):
    return char*length

def delim_string(length,string,char):
    len_part = int((length-len(string)-2) / 2)
    if len(char*len_part + ' ' + string + ' ' + char*len_part) < length:
        return char*len_part + ' ' + string + ' ' + char*len_part + char
    else:
        return char*len_part + ' ' + string + ' ' + char*len_part

##########################################################
################### Trajectory output ####################
##########################################################

def write_trajectory(aggs,BOX,r_ref,max_z,count_steps,k,write,modules,Material_selection):
    INDEX = 0
    particles = 0
    for agg in range(len(aggs)):
        if agg <= k:
            for PP in range(aggs[agg].Np):
                particles += 1
    write.write('ITEM: TIMESTEP\n%i\n' %(count_steps))
    write.write('ITEM: NUMBER OF ATOMS\n%i\n' %(particles))
    write.write('ITEM: BOX BOUNDS pp pp pp\n')
    write.write('%1.3e %1.3e\n' %(BOX[0] * r_ref, BOX[1] * r_ref))
    write.write('%1.3e %1.3e\n' %(BOX[2] * r_ref, BOX[3] * r_ref))
    write.write('0 %1.3e\n' %((max_z + 10) * r_ref))
    if 4 not in modules:
        write.write('ITEM: ATOMS id type x y z Radius\n')
        for agg in range(len(aggs)):
            if agg <= k:
                for PP in range(aggs[agg].Np):
                    write.write('%i %i %3.6f %3.6f %3.6f %3.6f\n' %(INDEX,agg, aggs[agg].x[PP] * r_ref,aggs[agg].y[PP] * r_ref,aggs[agg].z[PP] * r_ref, aggs[agg].r[PP] * r_ref))
                    INDEX += 1

    else:
        write.write('ITEM: ATOMS id type x y z Radius\n')
        for agg in range(len(aggs)):
            if agg <= k:
                for PP in range(aggs[agg].Np):
                    write.write('%i %i %3.6f %3.6f %3.6f %3.6f\n' %(INDEX,aggs[agg].TYPE, aggs[agg].x[PP] * r_ref,aggs[agg].y[PP] * r_ref,aggs[agg].z[PP] * r_ref, aggs[agg].r[PP] * r_ref))
                    INDEX += 1
