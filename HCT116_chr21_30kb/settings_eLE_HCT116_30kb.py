import numpy as np
#import inspect
import pandas as pd
import hoomd


##############################
# Simulation parameters
##############################
thermal_timesteps = 1e6
final_timesteps = 1.5e8
dump_period = 5e4
restart_period = 1e6


N_le = 10  # Number of loop extrusion agents/factors

def dump_fname(sim_id):
    return r'{sim_id}-dump.gsd'.format(sim_id=sim_id)
    
def restart_fname(sim_id):
    return r'{sim_id}-restart.gsd'.format(sim_id=sim_id)
    


# MD parameters
dt = 0.01
KT = 1

# Simulation box
N = 0  # do not change this number
N_pol = 930  
L0 = 55
L = np.array([L0, L0, L0])
Lhalf = 0.5*L
Dhalf = Lhalf
box = hoomd.data.boxdim(*L)

# LJ potential
sigma = 1.0
epsilon = 1.0
fene_r0 = 1.6
lj_rcut_rep = sigma*2**(1./6.)  # act on all particles from different types
# sigma*2.5
# act on all particles and their corresponding binders
lj_rcut_att = lj_rcut_rep
# average background lj potential
epsilon_background = 1.0  
lj_background_rcut = 1.5 

# loop-extrusion
LE_bond_type = 'harmonic'  # set to either 'harmonic' or 'fene'.
le_period = 5e2  # Number of MD steps to update/slide the loop extrusion factor
le_p_off = 0.0085  
le_extr_dist_thr = -1  # do not change it 
le_type = 'loop_extrusion'
le_K = 10
le_r0 = 1.1
# Number of MD steps to update/slide the loop extrusion factor
ctcf_df = pd.read_csv('HCT116_chr21_30kb/eLE_Peaks_HCT116.tsv', delim_whitespace=True)
stop_probs = ctcf_df[['forward', 'backward']].values
stop_probs[0, :] = 1.0
stop_probs[-1, :] = 1.0


##############################
# particle, bond and angle types
# read the particle information from a file
##############################
binder_bead_ratio = 0. 
typeid_fname = 'HCT116_chr21_30kb/polymer_HCT116.bed'
polymer_df = pd.read_csv(typeid_fname, delim_whitespace=True,
                         header=None, names=['pid', 'type'])
polymer_particle_types = list(set(polymer_df.type.values))
binder_particle_types = ['%s_binder' % x for x in polymer_particle_types]
particle_types = polymer_particle_types + binder_particle_types

bond_types = ['polymer', 'loop_extrusion']
le_typeid = bond_types.index(le_type)

snap_args = {'N': 0, 'box': box,
             'particle_types': particle_types, 'bond_types': bond_types}

# define polymers
polymers = []
pol_size = len(polymer_df)
for i in range(1):
    particles_typeid = np.array([particle_types.index(x)
                                 for x in polymer_df.type.values])
    #print ('particles_typeid : ', pol_size, particles_typeid)
    pol0 = dict(size=pol_size, bonds_typeid=0, particles_typeid=particles_typeid,
                pbc=np.array([True, True, True]), rcm_pos=np.zeros(3))
    polymers.append(pol0)

# define binders
binder_sizes = [int(binder_bead_ratio*np.sum(particles_typeid ==
                                             particle_types.index(x))) for x in polymer_particle_types]
#print('binder_sizes', binder_sizes)
binders = [dict(size=binder_sizes[binder_particle_types.index(
    bt)], particles_typeid=particle_types.index(bt)) for bt in binder_particle_types]
print(binders)
