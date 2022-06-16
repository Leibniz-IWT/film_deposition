""" Inputfile for the calculation of particle films
The names of the variables must not be changed! The order is irrelevant

##Abbrevations
PP = primary particle
agg = Aggregate
random = calculated by random value between two limits
"""
Version = '2.4.1'
#######################
# ----- General ----- #
#######################
# Names and formats used for all outputfiles
experiment = 'dda'
# Trajectory files
traj_format = '.trj'
img_format = '.pdf'
# Number of total particles
N_tot = 72

# Ballistic velocity
# Typical thermophoretic deposition velocities are 0.001-0,1 m/s-1
# Mädler et al. (2006) Sensors Actuators B 114, 283-95
Vel = 1e0
Mov_Method = 'Pe'          # Which parameter ist set? 'Pe' or 'Vel'
Pe = 1                      # Peclet-Number (Pe = r*v/D)
# Temperature in K
T = 298

# If Box_rel == True: length box = Box_size * r
# else length box = Box_size (nm)
Box_size = 280

# In which height above the highest particle will be the new aggregate initialized
# If Box_rel == True: it is also relative to r, else in nm
Height_above = 60

##########################
# ----- Parameters ----- #
##########################

# Distribution of N (none, random,lognorm, norm,predef))
DLA_Mode = 'CCAS'            # SA / CCAN (Clusters = number of clusters)/ CCAS (cluster = size of clusters)
N_dist = 'none'
Cluster = 6
N_min = 36                  # Minimum PP per agg (this is used if N_variable = False)
N_max = 0                   # Maximum PP per agg
N_sigma = 0.233             # Standard derivation in case of norm or lognorm
# Distribution of PP radius (none, random, lognorm, norm, predef)
# Consider that predef needs rmin to make it dimensionless and see predef_file
r_dist = 'lognorm'
rmin = 4.479                # Minimum PP diameter in nm
#(is used as mean for none, lognorm and norm)
rmax = 7.5                  # Maximum PP diameter in nm (only for random)
r_sigma = 0.370             # Standard derivation of distribution (lognorm, norm)

# Typical values here are D_f = 1.8 and k_f = 1.3 (Friedlander 2000 'Smoke, Dust and Haze', p. 227)
# or Df = 1.78 +- 0.1/ kf = 1.3+-0.2 (Sorensen 2011 Aerosol Science and Technoloy, 45:765-779)
D_f = 1.8                   # Fractal dimension of the aggs
k_f = 1.3                   # Prefactor value for fractal dimension
epsilon = 1.001             # Factor for minimum distance between two PP in one agg
delta = 1.01                # Factor for maximum distance between two PP in one agg
deformation = 1.00          # Artificial deformation of the particle

dt = 1                      # The timestep is chosen due to the RMSD. It can be adjusted here

analysis_height = 0.9       # This relative height is used to evaluate porosity and percolation
# so that single Aggregates reaching out are omitted

# File containing the diameter and propability of the predefined size distribution:
# rp (in nm) propability effective dp, eg:
# 1.5     0.12    3.06
# 2.0     0.24    4.08
# ...

# The files for manual distributions have to be in the DLA_Source folder
pre_dp_file = 'table.polydisperse_FSP_SF_Medium'
use_effective_dp = True     # With this in the granular file the effective diameter defined in
# pre_dp_file is used instead of actual diameter
pre_N_file = 'table.lognorm_N_dist_CPS'
######################
# ----- Modules----- #
######################
''' Here you can choose from several options, which change the overall behavior '''
modules = [2, 3]
'''
1. Multiple deposition (Aggregates rotates after ist first deposition until it has hit 3 times)
2. Calculate Coordination
3. Calculate Percolation
4. Multicomponent system

'''
''' How will be chosen which material an aggregate consists of:
    0 Always rho_1
    1 Change for each aggregate
    2 Random distribution
'''
mixingratio = 0.5
# Threshold sets the max distance until a particle
# is count as neighbor
'''
Threshold_mode:
1 One rp distance = neighbor (polydisperse: mean r_0,r_1)
2 delta (as above) will be multiplied
3 custom value in nm (threshold)
'''
coord_threshold_mode = 2
coord_threshold = 1
perc_threshold_mode = 1
perc_threshold = 6

Material_selection = 1
# Density of the particles in g/m3: TiO2: 4.23e3 // SnO2: 6.85e3
''' How will be chosen which material an aggregate consists of:
    0 Always rho_1
    1 Change for each aggregate
    2 Random distribution
'''
rho_1 = 4.23e3
rho_2 = 6.85e3
#####################
# ----- Output----- #
#####################
""" Here you can choose which outputfiles will be created.
multiple choices are possible

These two file will always be saved:
A) .txt file with all the coordinates (in m) (index,agg,r,x,y,z) + a .txt with all aggregate and simulation data (velocity,box,D_f,Pe,rho,part.-dist type, diff.coeff, mass, PP, etc.)
B) trajectory file with the final geometry of the particle film for e g. ovito
"""
output = [1, 2, 3, 5, 7, 8]
'''
These files are optional:
1 log file
2 LIGGGHTS data file hybrid granular/molecular (head,index,type,diameter,density,x,y,z, molecule(aggregate))
3 LIGGGHTS data file granular (head,index,type,diameter,density,x,y,z)
4 Save the whole trajectory of the deposition
5 VTK outputfile for paraview
6 trajectory file with the format for mixing_states
7 A histogram of the primary particle size
8 A histogram of the amount of primary particles per Aggregate
9 Building Time of kd-tree (neighbor list)
10 LAMMPS data file molecular (head,index,agg,x,y,z,0,0,0)
'''
# After how many aggregates an output in the log file is written
log_step = 10
# Set the stepsize for '4'
steps_trajectory = 1

#######################
# ----- Methods ----- #
#######################
# Friction Coefficient:     1 = Chen and Dahneke, 2 = Cunningham, 3 = Zhang (read from prepared file diffusionfile)
Fric_method = 1
diffusionfile = '/home/baric/src/git/dda/table.diffusion_zhang'
# Use method to calculate geometrical instead of mass center of aggregate
geom_center = 0
# Calculation of mean radius. 0 = rmin, 1 = arithmetic, 2 = geometric, 3 = harmonic, 4 = sauter
# Acc. to Eggersdorfer et al. 2012 Aerosol Science and Technology 46:3, 347-353: mean_r = 2
# Acc. to Eggersdorfer et al. 2012 Journal of Colloid and Interface Science 387: 12-23, mean_r = 4
mean_r = 2
# Method to calculate a random standard derivation. 1 = Box,Muller 1958, 2 = Mädler 2006
sdev_method = 2

# How often will it be tried to set a single PP within an agg
iteration_limit = 1e4

# The angle increment that is used for each rotation step
'''
Warning! A high value results in a very high resolution trajectory. Consider this!
'''
ang_inc = 72
# How often can the direction of a rotation change before it stops
rot_attemps = 100

# Folder used to store outputdata
output_folder = './'


# Folder on Ranjid: '/home/baric/Piko/DLA/DLA_source/input_archive/'
# Folder on Baric-PC: '/home/baric/Dokumente/Simulationen/DLA/Archive_Inputfiles/'
# Folder on Mac: '/home/valentinbaric/Tools/DLA/Archiv_Input/'

########################
#### Misc ##############
########################
send_email = True
toaddr = 'v.baric@iwt.uni-bremen.de'
extra_text = ''

########################
#### Debugging Mode ####
########################

# Save the last particle as a new standard particle
create_standard = False
debugging = [2]
""" Here you can activate several debugging options:
!!! If you want no debugging, you have to insert 0 !!!
    The order is not relevant, you can enter the numbers you want

    1 Only the generation of the aggregates, no transport
    2 start at defined coordinate (see below)
    3 start in an area defined by 'define_area'
    4 search for overlapping particles including boundaries (after deposition)
    5 Check overlap after every aggregate
    6 stop simulation in case of overlap (together with (13))
    7 use the same seed for random number generator (custom_seed)
    8 use a custom seed for random number generator aggregates
    9 suppress progress bar
"""
# In case of '10' this is the starting coordinate (consider that the maximum coordinate is BOX_size/2)
debug_coord = [0, 0]
define_area = 100
custom_seed_traj = 77777
custom_seed_aggs = 88888
boundaries = False
# If a custom seed for aggregates is selected, all aggregates will look be same
# if a monodisperse number distribution is selected!
