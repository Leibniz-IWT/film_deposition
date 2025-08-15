import os, glob, csv, math
import matplotlib.pyplot as plt
import polyscope as ps
import numpy as np
from classes.particlelayer import *
import numpy as np
import json
import polyscope as ps
from scipy.optimize import shgo


# Data imports:
file  = './case_studies/hetero/results/hetero.trj'
#file  = './case_studies/hetero/resultshetero.trj'
layer = particlelayer(file, nodescription=False).layer

# Shape
layer._data.shape

# Extract positions in space
points = np.zeros([layer.data.shape[0], 3])
points[:, 0] = layer.data[:, layer.lib['x']]
points[:, 1] = layer.data[:, layer.lib['y']]
points[:, 2] = layer.data[:, layer.lib['z']]
particles = points  # used later
# Extract the radius
radii  = np.zeros(layer.data.shape[0])
radii = layer.data[:, layer.lib['r']]

if 0:
    # Plot in polyscope
    ps.init()
    ps.set_up_dir("z_up")
    ps_cloud = ps.register_point_cloud("Particles", points)
    ps_cloud.add_scalar_quantity("Particle radii", radii)
    ps_cloud.set_point_radius_quantity("Particle radii", autoscale=False)
    ps_cloud.set_color((1.0, 1.0, 1.0))
    ps.show()


# Util functions
def compute_contact_stats(
        pos_A, pos_B,                # (NA,3)  & (NB,3)
        rad_A, rad_B,                # (NA,)   & (NB,)
        tol      = 0.01,             # fractional slack on touch distance
        block_sz = 4096              # tune if memory is tight
    ):
    """
    Returns:
        homo_contacts, hetero_contacts, coord_homo, coord_hetero
    (See earlier doc-string for full description.)
    """
    # ------------------------------------------------------------------
    # 0.  Merge arrays for easier indexing
    # ------------------------------------------------------------------
    pos   = np.vstack([pos_A, pos_B])
    radii = np.hstack([rad_A, rad_B])
    ptype = np.hstack([np.zeros(len(pos_A), dtype=np.int8),
                       np.ones (len(pos_B), dtype=np.int8)])
    N     = len(radii)

    homo_contacts   = 0
    hetero_contacts = 0
    coord_homo      = np.zeros(N, dtype=np.uint32)
    coord_hetero    = np.zeros(N, dtype=np.uint32)

    # ------------------------------------------------------------------
    # 1.  Pair search – inter-block **and** intra-block
    # ------------------------------------------------------------------
    for lo in range(0, N, block_sz):
        hi = min(lo + block_sz, N)

        # --------------------------------------------------------------
        # 1a.  Pairs where j ≥ hi   (different blocks)   ← as before
        # --------------------------------------------------------------
        if hi < N:                                        # skip if last block
            P_block   = pos[lo:hi]
            radii_blk = radii[lo:hi][:, None]

            diff   = pos[hi:] - P_block[:, None, :]
            dist2  = np.sum(diff**2, axis=-1)
            sum_r  = radii_blk + radii[hi:]
            thresh2 = (sum_r * (1.0 + tol))**2

            ii, jj = np.where(dist2 <= thresh2)
            i_idx  = lo + ii
            j_idx  = hi + jj

            same = (ptype[i_idx] == ptype[j_idx])
            homo_contacts   += int(np.count_nonzero(same))
            hetero_contacts += int(np.count_nonzero(~same))

            for a, b in zip(i_idx, j_idx):
                if same[np.where(i_idx==a)[0][0]]:
                    coord_homo[a]   += 1
                    coord_homo[b]   += 1
                else:
                    coord_hetero[a] += 1
                    coord_hetero[b] += 1

        # --------------------------------------------------------------
        # 1b.  Pairs *inside* [lo:hi)  (same block) – upper triangle
        # --------------------------------------------------------------
        if hi - lo > 1:                                   # need ≥2 points
            P_block   = pos[lo:hi]
            r_block   = radii[lo:hi]

            # full pairwise distances inside the block
            diff      = P_block[:, None, :] - P_block[None, :, :]
            dist2     = np.sum(diff**2, axis=-1)

            # radii sums matrix
            sum_r     = r_block[:, None] + r_block[None, :]
            thresh2   = (sum_r * (1.0 + tol))**2

            # keep only upper-triangle (k=1 removes diagonal)
            iu, ju = np.triu_indices(hi - lo, k=1)
            mask   = dist2[iu, ju] <= thresh2[iu, ju]
            if not np.any(mask):
                continue

            i_idx = lo + iu[mask]
            j_idx = lo + ju[mask]

            same = (ptype[i_idx] == ptype[j_idx])
            homo_contacts   += int(np.count_nonzero(same))
            hetero_contacts += int(np.count_nonzero(~same))

            for a, b in zip(i_idx, j_idx):
                if ptype[a] == ptype[b]:
                    coord_homo[a]   += 1
                    coord_homo[b]   += 1
                else:
                    coord_hetero[a] += 1
                    coord_hetero[b] += 1

    return homo_contacts, hetero_contacts, coord_homo, coord_hetero


base_dirs = [
    "./case_studies/hetero/Streamline_PBE/DFSP_new_left",
    "./case_studies/hetero/Streamline_PBE/DFSP_new_right",
]

streamlines = {}          # {global_idx: {col: [values,…]}}
global_idx  = 0           # will increment for every file we load

for base in base_dirs:
    pattern = os.path.join(base, "processed_streamline_*.csv")
    for fp in sorted(glob.glob(pattern)):          # keep a stable order
        with open(fp, newline='') as f:
            reader = csv.DictReader(f)
            data   = {h: [] for h in reader.fieldnames}
            for row in reader:
                for h in reader.fieldnames:
                    try:
                        data[h].append(float(row[h]))
                    except (ValueError, TypeError):
                        data[h].append(math.nan)

        streamlines[global_idx] = data
        global_idx += 1

print(f"Loaded {len(streamlines)} streamlines total")
###
# plotting for validation
plot_keys = ["O2A", "O2B", "T"]
time_key  = "IntegrationTime"

for key in plot_keys:
    #plt.figure()
    for idx in sorted(streamlines):
        d = streamlines[idx]
        #plt.loglog(d[time_key], d[key], alpha=0.25, linewidth=0.8)
    #plt.xlabel("Real time (s)")
    #plt.ylabel(key)
    #plt.title(f"{key} across {len(streamlines)} streamlines (Left + Right)")
    #plt.tight_layout()

#plt.show()
###
# Quick tests
T = d['T']  # d is the last data set we used in the loop
t_n_start_i = np.argmax(T)
t_n_start_i = 0  # New code use all data
t_n_start_i, len(T)


# Kinetics

# Analytical solutions of mean coordination numbers:
homo_cn_anal = 4.0
#>>> 4.0**(2/3)  # 2D --> 3D
#2.5198420997897464
homo_cn_anal = 2.5198420997897464
homo_cn_anal = 2.33
homo_cn_anal = 2.33**0.5

heter_cn_anal = 0.4
#>>> 0.4**(2/3)  # 2D --> 3D
#0.5428835233189814
heter_cn_anal = 0.5428835233189814
heter_cn_anal = 0.7
heter_cn_anal = 0.7**0.5



# Define KMC rate laws and parameters
# Parameters (to be optimized)
k_A = 1  #  Arrhenius prefactor for compound A
k_B = 1  #  Arrhenius prefactor for compound B
Ea_A = 1  # Activation energy for compound A (should be based on nucleation)
Ea_B = 1  # Activation energy for compound B (should be based on nucleation)

# Parameters (fixed)
n_A = 1  # Rate law order of compound A
n_B = 1  # Rate law order of compound B
R = 8.314  # Universal gas constant J kmol-1 K-1

# Rate law
def r_i(c, T, k, E):
    return k * np.exp(-E / (R * T)) * c

obj = 0.0
print(f'Running quick test of data for single optimization call...')
for idx in sorted(streamlines):
    # Time discretisation (here we determine the real time spent available for nucleation)
    d = streamlines[idx]
    tn = np.array(d['IntegrationTime'][t_n_start_i:]) - np.array(d['IntegrationTime'][t_n_start_i])  # time - starting time
    Tn = np.array(d['T'][t_n_start_i:])
    CA_n = np.array(d['O2A'][t_n_start_i:])
    CB_n = np.array(d['O2B'][t_n_start_i:])
    #if Tn.shape[0] > CA_n.shape[0]:
    #    CA_n = np.array(d['O2A'][t_n_start_i:Tn.shape[0]])
    #elif CA_n.shape[0] < Tn.shape[0]:
    #    CA_n = np.array(d['O2A'][t_n_start_i:Tn.shape[0]])
    #CB_n = np.array(d['O2B'][t_n_start_i:Tn.shape[0]])

    # Usable values (print in debugging:
    #np.sum(tn), len(tn), len(Tn), len(CA_n), len(CB_n)
    #particles.shape[0]/len(tn)

    # Assume even split for now (nucleation rates probably based on cooling rate, which is relatively linear?
    p_formed = np.zeros_like(tn)  # Number of particles nucleated at every time step
    p_formed[:] = int(particles.shape[0]/len(tn)) # Can be modified here later for more realistic distributions.
    p_formed[-1] =  particles.shape[0] - np.sum(p_formed[:-1])
    #p_formed, np.sum(p_formed)

    # Reaction rates
    r_A_range = r_i(CA_n, Tn, k_A, Ea_A)
    r_B_range = r_i(CB_n, Tn, k_B, Ea_A)
    r_t = r_A_range + r_B_range
    # KMC relative formation rates
    ra = r_A_range / r_t
    rb = r_B_range / r_t
    #ra, rb

    # Main nucleation loop
    # The logic is that each time step we used the KMC rates derived above to
    # find the nucleating particles each time step
    particles_A_i = []  # Indices of particle A
    particles_B_i = []  # Indices of particle B
    # KMC formation each time step
    particles_i = 0  # current particle index tracker
    for i, pf in enumerate(p_formed):
        for pi in range(int(pf)):
            # Distributed in kmc ration
            if np.random.rand() < ra[i]:
                particles_A_i.append(particles_i)
            else:
                particles_B_i.append(particles_i)

            particles_i += 1  # Current particle index

    # Post processing:
    Particles_A = particles[particles_A_i]
    Particles_B = particles[particles_B_i]
    radii_A = radii[particles_A_i]
    radii_B = radii[particles_B_i]
    # Compute contacts
    hc, hetc, c_homo, c_hetero = compute_contact_stats(
        Particles_A, Particles_B, radii_A, radii_B)   # try 10 % slack
    # Compute mean coordination numbers
    homo_cn = hc / (Particles_A.shape[0] + Particles_B.shape[0])
    heter_cn = hetc / (Particles_A.shape[0] + Particles_B.shape[0])
    #print(f'homo_cn = {homo_cn}')
    #print(f'heter_cn = {heter_cn}')
    # Compute obj by norm to experimental values
    obj += np.linalg.norm(homo_cn - homo_cn_anal) + 1e3 * np.linalg.norm(heter_cn - heter_cn_anal)

len(particles_A_i), len(particles_B_i)
print(f'Done: obj = {obj}')
#print(f'obj = {obj}')

import numpy as np

# Assumed globals (you must ensure these are defined before optimization)
# - streamlines: dict of streamlines with keys 'IntegrationTime', 'T', 'O2A', 'O2B'
# - particles: ndarray of particle positions
# - radii: ndarray of particle radii
# - t_n_start_i: int, start index for nucleation time
# - compute_contact_stats: function(Particles_A, Particles_B, radii_A, radii_B)

def objective(params, homo_cn_anal=0.7, heter_cn_anal=0.25, n_A=1,  n_B=1):
    # Unpack optimization parameters
    k_A, k_B, Ea_A, Ea_B = params
    #k_A, k_B, Ea_A, Ea_B, n_A, n_B  = params

    # Fixed parameters (note, not actually added to equations yet)
    n_A = n_A
    n_B = n_B
    R = 8.314  # Universal gas constant

    # Analytical target coordination numbers
    homo_cn_anal = homo_cn_anal # homo_cn_anal 2D = 4.0
    heter_cn_anal = heter_cn_anal # heter_cn_anal = 0.4

    # Rate law function
    def r_i(c, T, k, E):
        return k * np.exp(-E / (R * T)) * c

    # Objective value
    obj = 0.0

    for idx in sorted(streamlines):
        d = streamlines[idx]

        tn = np.array(d['IntegrationTime'][t_n_start_i:]) - d['IntegrationTime'][t_n_start_i]
        Tn = np.array(d['T'][t_n_start_i:])
        CA_n = np.array(d['O2A'][t_n_start_i:])
        CB_n = np.array(d['O2B'][t_n_start_i:])

        # Nucleation event distribution
        p_formed = np.zeros_like(tn, dtype=int)
        p_formed[:] = int(particles.shape[0] / len(tn))
        p_formed[-1] = particles.shape[0] - np.sum(p_formed[:-1])

        # Reaction rates
        r_A_range = r_i(CA_n, Tn, k_A, Ea_A)
        r_B_range = r_i(CB_n, Tn, k_B, Ea_B)
        r_t = r_A_range + r_B_range

        # Avoid divide-by-zero
        r_t[r_t == 0] = 1e-12
        ra = r_A_range / r_t
        rb = r_B_range / r_t

        # KMC particle assignment
        particles_A_i = []
        particles_B_i = []
        particles_i = 0
        for i, pf in enumerate(p_formed):
            for _ in range(int(pf)):
                if np.random.rand() < ra[i]:
                    particles_A_i.append(particles_i)
                else:
                    particles_B_i.append(particles_i)
                particles_i += 1

        Particles_A = particles[particles_A_i]
        Particles_B = particles[particles_B_i]
        radii_A = radii[particles_A_i]
        radii_B = radii[particles_B_i]

        hc, hetc, _, _ = compute_contact_stats(Particles_A, Particles_B, radii_A, radii_B)

        total_particles = Particles_A.shape[0] + Particles_B.shape[0]
        homo_cn = hc / total_particles if total_particles > 0 else 0
        heter_cn = hetc / total_particles if total_particles > 0 else 0
        #print(f'homo_cn = {homo_cn }')
        #print(f' homo_cn_anal = { homo_cn_anal}')
        obj += np.linalg.norm(homo_cn - homo_cn_anal) + np.linalg.norm(heter_cn - heter_cn_anal)

    return obj



#%%
# Run the optimization
#%%

if 0:
    # Define bounds for each parameter: (k_A, k_B, Ea_A, Ea_B)
    bounds = [
        (1e-3, 100),  # k_A
        (1e-3, 100),  # k_B
        (1e-1, 1000), # Ea_A
        (1e-1, 1000)  # Ea_B
    ]

    #args=(homo_cn_anal=0.7, heter_cn_anal=0.25, n_A=1,  n_B=1)
    result = shgo(objective, bounds,
                  #n=1000
                  )
    print(result)

# Define bounds for each parameter: (k_A, k_B, Ea_A, Ea_B)
bounds = [
    #(1e-3, 10000),   # k_A
    #(1e-5, 10000),   # k_A
    #(1e-5, 10),   # k_A
    (1e-9, 1e-5),   # k_A
    #(1e-3, 10000),   # k_B
    #(1e-5, 10000),   # k_B
    #(1e-5, 10),   # k_B
    (1e-9, 1e-5),   # k_B
    (1e0, 1e4),  # (1e-2, 1e8),  # Ea_A  # 10^7    Literature ~ 1.8 x 10^5 J mol-1 K-1  --> 10^2 J kmol-1 K-1
    (1e0, 1e4)  # (1e-2, 1e8)   # Ea_B
]

n=5000
#sampling_method='sobol'
sampling_method='simplicial'
#args=(homo_cn_anal=0.7, heter_cn_anal=0.25, n_A=1,  n_B=1)
args=(homo_cn_anal, heter_cn_anal, 1, 1)
#options = {'maxiter': 15000}
minimizer_kwargs = {'method':'Nelder-Mead'}
print('Starting optimization run...')
print(f'Running with n={n} and bounds = {bounds}')
print(f'sampling_method = {sampling_method}')
result = shgo(objective, bounds,
              #n=1000,
              n=n,
              sampling_method=sampling_method,
              #options=options,
              minimizer_kwargs=minimizer_kwargs,
              workers=32,  # ~4-10 hours+
              args=args
              )



print(result)
#%%

"""
Results

 python hetero_kinetics_DFSP.py 

### Reading Particle Structure File
# Filename: ./case_studies/hetero/results/hetero.trj
    File ./case_studies/hetero/results/hetero.trj using units nano and type dda
# File Format: trj
# Units: nano
# Primary Particles: 2,000
# Aggregates: 19
# Box Width: 1.500e+02 x 1.500e+02
# Film Height: 2.974e+02
### Done


Loaded 138 streamlines total
Running quick test of data for single optimization call...
Done: obj = 189.5559414718399
Starting optimization run...
Running with n=5000 and bounds = [(1e-08, 1e-05), (1e-08, 1e-05), (10.0, 1000.0), (10.0, 1000.0)]
sampling_method = simplicial
 message: Optimization terminated successfully.
 success: True
     fun: 189.55594147183987
    funl: [ 1.896e+02  1.896e+02 ...  1.896e+02  1.896e+02]
       x: [ 6.344e-07  9.376e-06  1.956e+02  3.194e+02]
      xl: [[ 6.344e-07  9.376e-06  1.956e+02  3.194e+02]
           [ 4.600e-06  9.376e-06  3.194e+02  5.669e+02]
           ...
           [ 6.344e-07  6.921e-06  1.956e+02  4.431e+02]
           [ 1.002e-05  1.002e-05  1.002e+03  9.969e+02]]
     nit: 1
    nfev: 14266
   nlfev: 9264
   nljev: 0
   nlhev: 0


$ python hetero_kinetics_DFSP.py 

### Reading Particle Structure File
# Filename: ./case_studies/hetero/results/hetero.trj
    File ./case_studies/hetero/results/hetero.trj using units nano and type dda
# File Format: trj
# Units: nano
# Primary Particles: 2,000
# Aggregates: 19
# Box Width: 1.500e+02 x 1.500e+02
# Film Height: 2.974e+02
### Done


Loaded 138 streamlines total
Running quick test of data for single optimization call...
Done: obj = 281.5890000000005
Starting optimization run...
Running with n=500 and bounds = [(1e-08, 1e-05), (1e-08, 1e-05), (10.0, 1000.0), (10.0, 1000.0)]
sampling_method = simplicial
 message: Failed to find a feasible minimizer point. Lowest sampling point = 281.5890000000005
 success: False
     fun: 281.5890000000005
       x: [ 1.000e-08  1.000e-08  1.000e+01  1.000e+01]
     nit: 1
    nfev: 500
   nlfev: 0
   nljev: 0
   nlhev: 0
    tnev: 500



$ python hetero_kinetics_DFSP.py 

### Reading Particle Structure File
# Filename: ./case_studies/hetero/results/hetero.trj
    File ./case_studies/hetero/results/hetero.trj using units nano and type dda
# File Format: trj
# Units: nano
# Primary Particles: 2,000
# Aggregates: 19
# Box Width: 1.500e+02 x 1.500e+02
# Film Height: 2.974e+02
### Done


Loaded 138 streamlines total
Running quick test of data for single optimization call...
Done: obj = 189.55594147183987
Starting optimization run...
Running with n=5000 and bounds = [(1e-08, 1e-05), (1e-08, 1e-05), (10.0, 1000.0), (10.0, 1000.0)]
sampling_method = simplicial
^[[A^[[A^[[DMNerg
 message: Optimization terminated successfully.
 success: True
     fun: 189.55594147183987
    funl: [ 1.896e+02  1.896e+02 ...  1.896e+02  1.896e+02]
       x: [ 2.508e-06  1.000e-05  1.000e+01  2.576e+02]
      xl: [[ 2.508e-06  1.000e-05  1.000e+01  2.576e+02]
           [ 1.000e-08  1.000e-05  7.526e+02  1.001e+03]
           ...
           [ 5.068e-06  1.274e-06  2.607e+02  4.924e+02]
           [ 8.330e-06  6.344e-07  3.194e+02  5.669e+02]]
     nit: 1
    nfev: 12325
   nlfev: 7323
   nljev: 0
   nlhev: 0


$ python hetero_kinetics_DFSP.py 

### Reading Particle Structure File
# Filename: ./case_studies/hetero/results/hetero.trj
    File ./case_studies/hetero/results/hetero.trj using units nano and type dda
# File Format: trj
# Units: nano
# Primary Particles: 2,000
# Aggregates: 19
# Box Width: 1.500e+02 x 1.500e+02
# Film Height: 2.974e+02
### Done


Loaded 138 streamlines total
Running quick test of data for single optimization call...
Done: obj = 83755.98951951254
Starting optimization run...
Running with n=5000 and bounds = [(1e-09, 1e-05), (1e-09, 1e-05), (1.0, 10000.0), (1.0, 10000.0)]
sampling_method = simplicial
 message: Optimization terminated successfully.
 success: True
     fun: 189.55594147183987
    funl: [ 1.896e+02  1.896e+02 ...  1.896e+02  1.896e+02]
       x: [ 6.563e-06  8.438e-06  5.938e+03  3.524e+03]
      xl: [[ 6.563e-06  8.438e-06  5.938e+03  3.524e+03]
           [ 9.375e-06  6.875e-06  3.126e+03  1.876e+03]
           ...
           [ 8.860e-06  5.938e-06  3.438e+03  3.438e+03]
           [ 2.509e-06  7.524e-06  5.016e+03  7.594e+03]]
     nit: 1
    nfev: 19200
   nlfev: 14200
   nljev: 0
   nlhev: 0
"""
