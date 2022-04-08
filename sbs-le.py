#! /usr/bin/env python3

import hoomd
import numpy as np
import itertools
import os
import sys
from hoomd import md
import argparse
import sbs_class
import saw
import le4hoomd

parser = argparse.ArgumentParser()
parser.add_argument('--settings',
                    type=str,
                    help='the setting file',
                    default='default-settings.py',
                    dest='settings')
parser.add_argument('--sim-id',
                    type=str,
                    help='the simulaiton id',
                    default='',
                    dest='sim_id')
parser.add_argument('--resume',
                    type=bool,
                    help='resume (default: False)',
                    default=False,
                    dest='resume')
parser.add_argument('--skip-rm-overlaps',
                    type=bool,
                    help='skip rmoving overlaps (default: False)',
                    default=False,
                    dest='skip_rm_overlaps')
parser.add_argument('--skip-thermalization',
                    type=bool,
                    help='skip thermalization (default: False)',
                    default=False,
                    dest='skip_thermalization')
parser.add_argument(
    '--init-file',
    type=str,
    help=
    'the initial configuration file. If pass nothing, the restart file will be checked.',
    default='init-restart.gsd',
    dest='init_file')
parser.add_argument(
    '--Nle',
    type=int,
    default=-1,
    help='Number of LEF. If nothing passed, settings.N_le is used.')
parser.add_argument('--hoomd',
                    '--hoomd_args',
                    default='gpu',
                    dest='hoomd_args',
                    type=str)
args = parser.parse_args()
print(args.settings)
exec ('import {} as settings'.format(args.settings[:-3]))

######################
## Read the polymer ##
## bead/binder      ##
######################
polymer_particle_types = settings.polymer_particle_types
binder_particle_types = settings.binder_particle_types
particle_types = settings.particle_types

if args.Nle == -1:
    N_le = settings.N_le
else:
    N_le = args.Nle
########################
## Simulation context ##
########################
sbs = sbs_class.sbs(hoomd, settings)
sbs.init_world(mode=args.hoomd_args)

dump_fname = settings.dump_fname(args.sim_id)
restart_fname = settings.restart_fname(args.sim_id)

restart = args.resume  # flag the simulation mode!
if (args.init_file != 'init-restart.gsd') or args.resume:
    # read the restart file
    if restart:
        system = hoomd.init.read_gsd(filename=restart_fname,
                                     restart=args.init_file)
    else:
        system = hoomd.init.read_gsd(filename=args.init_file, time_step=0)
else:
    snapshot = hoomd.data.make_snapshot(**settings.snap_args)
    for pol0 in settings.polymers:
        print(pol0)
        sbs.add_random_polymer(snapshot,
                               L=settings.L * np.array([1, 1, 0.9]),
                               **pol0)
    for binders0 in settings.binders:
        if binders0['particles_typeid'] == particle_types.index('A_binder'):
            continue
        sbs.add_random_particles(snapshot, L=settings.L * 0.95, **binders0)
    # read the first snapshot
    system = hoomd.init.read_snapshot(snapshot)

# Neighborlist
nl = md.nlist.cell(check_period=1)
nl.set_params(r_buff=1.1)

#####################
## pairs and bonds ##
#####################
sigma = settings.sigma
lj_rcut_rep = settings.lj_rcut_rep
lj_rcut_att = settings.lj_rcut_att
lj_background_rcut = settings.lj_background_rcut
fene_r0 = settings.fene_r0

dpd = hoomd.md.pair.dpd(r_cut=sigma, kT=0.0, nlist=nl, seed=1)
lj_rep = hoomd.md.pair.lj(r_cut=lj_rcut_rep, nlist=nl)
lj_att = hoomd.md.pair.lj(r_cut=lj_rcut_att, nlist=nl)
lj_background = hoomd.md.pair.lj(r_cut=lj_background_rcut, nlist=nl)
lj_rep.set_params(mode='shift')
lj_att.set_params(mode='shift')
lj_background.set_params(mode='shift')

bond_fene = hoomd.md.bond.fene()
bond_fene.bond_coeff.set(type='polymer',
                         k=30.0 * 1,
                         r0=fene_r0 * sigma,
                         sigma=sigma,
                         epsilon=1)
bond_fene.bond_coeff.set(type='loop_extrusion',
                         k=0 * 1,
                         r0=fene_r0 * sigma,
                         sigma=sigma,
                         epsilon=1)
if settings.LE_bond_type == 'harmonic':
    bond_le = hoomd.md.bond.harmonic()
    bond_le.bond_coeff.set(type='polymer', k=0.0, r0=1.6 * sigma)
    bond_le.bond_coeff.set(type='loop_extrusion',
                           k=settings.le_K,
                           r0=settings.le_r0 * sigma)
elif settings.LE_bond_type == 'fene':
    bond_le = hoomd.md.bond.fene()
    bond_le.bond_coeff.set(type='polymer',
                           k=0.0,
                           r0=2.0 * sigma,
                           sigma=sigma,
                           epsilon=1),
    bond_le.bond_coeff.set(type='loop_extrusion',
                           k=settings.le_K,
                           r0=settings.le_r0 * sigma,
                           sigma=sigma * 1.,
                           epsilon=1)

####################
## init gsd files ##
####################
group_all = hoomd.group.all()
group_polymer = hoomd.group.tags(name='polymer',
                                 tag_min=0,
                                 tag_max=settings.pol_size - 1)
dump_gsd = hoomd.dump.gsd(
    dump_fname,
    period=settings.dump_period,
    group=group_all,
    overwrite=False,
    phase=0,
    dynamic=['attribute', 'property', 'momentum', 'topology'])
restart_gsd = hoomd.dump.gsd(
    restart_fname,
    period=settings.restart_period,
    group=group_all,
    overwrite=True,
    truncate=True,
    dynamic=['attribute', 'property', 'momentum', 'topology'])
####################
## remove_overlap ##
####################
dpd.pair_coeff.set(particle_types, particle_types, A=30, gamma=2.0)
lj_rep.disable()
lj_att.disable()
lj_background.disable()

if hoomd.get_step() >= settings.final_timesteps:
    sys.exit()
if not (args.skip_rm_overlaps):  # Check if it is really needed TODO
    sbs.remove_overlaps(group=group_all)
dpd.disable()

####################
## thermalization ##
####################
lj_rep.enable()
lj_rep.pair_coeff.set(particle_types,
                      particle_types,
                      epsilon=1,
                      sigma=settings.sigma,
                      r_cut=settings.lj_rcut_rep)
normal_md = md.integrate.mode_standard(dt=settings.dt)
langevin = md.integrate.langevin(group=group_all,
                                 kT=settings.KT,
                                 seed=np.random.randint(1, 100000 + 1))
for p0 in particle_types:
    langevin.set_gamma(p0, gamma=2.0)

print('restart:', restart, 'resume:', args.resume)
if not (restart):
    if not (args.skip_thermalization):
        hoomd.run_upto(int(settings.thermal_timesteps * 1))
    # loop-extrusion agents
    if N_le > 0:
        le_agents = le4hoomd.LEagents(system=system,
                                      stop_probs=settings.stop_probs,
                                      N_le=N_le,
                                      p_off=settings.le_p_off,
                                      move_dist_thr=settings.le_extr_dist_thr,
                                      update_period=settings.le_period,
                                      CTCF_type='random')
        le_agents.le_type = 'loop_extrusion'
        le_agents.le_typeid = settings.bond_types.index(le_agents.le_type)
        loop_extrusion = hoomd.analyze.callback(
            callback=le_agents.update_loop_extrusion2,
            period=settings.le_period)
        snap = system.take_snapshot(all=True)
        ctcf = le_agents.stop_probs.copy()
        ctcf[:, 1] *= 1.5  # just for visualization
        snap.particles.charge[:settings.pol_size] = np.sum(ctcf, axis=1)
        snap.bonds.resize(snap.bonds.N + N_le)
        le_indices = np.random.choice(settings.pol_size - 1,
                                      N_le,
                                      replace=False)
        for i, a1 in enumerate(le_indices):
            snap.bonds.group[-i - 1] = [a1, a1 + 1]
            snap.bonds.typeid[-i - 1] = settings.le_typeid
            print(i + 1, a1, a1 + 1)
        print('bonds = %d' % snap.bonds.N)
        np.savetxt(dump_fname[:-4] + '-ctcf.txt',
                   le_agents.stop_probs,
                   fmt='%.2f',
                   header='backward \t forward')
        system.restore_snapshot(snap)

if restart:
    # loop-extrusion agents
    if N_le > 0:
        ctcf = np.genfromtxt(dump_fname[:-4] + '-ctcf.txt')
        le_agents = le4hoomd.LEagents(system=system,
                                      stop_probs=ctcf,
                                      N_le=N_le,
                                      p_off=settings.le_p_off,
                                      move_dist_thr=settings.le_extr_dist_thr,
                                      update_period=settings.le_period,
                                      CTCF_type='all')
        le_agents.le_type = 'loop_extrusion'
        le_agents.le_typeid = settings.bond_types.index(le_agents.le_type)
        loop_extrusion = hoomd.analyze.callback(
            callback=le_agents.update_loop_extrusion2,
            period=settings.le_period)

##############################
## normal MD: attraction ##
##############################
lj_background.enable()
lj_att.enable()
# First remove interaction between all pairs of particles
lj_background.pair_coeff.set(particle_types,
                             particle_types,
                             sigma=1.0,
                             epsilon=settings.epsilon_background,
                             r_cut=False)
lj_background.pair_coeff.set(
    polymer_particle_types,
    binder_particle_types,
    sigma=1.0,
    epsilon=settings.epsilon_background,
    r_cut=lj_background_rcut)  # Then, just define it between beads and binders

lj_att.pair_coeff.set(particle_types,
                      particle_types,
                      epsilon=settings.epsilon,
                      sigma=1.0,
                      r_cut=False)
for p in polymer_particle_types:
    lj_background.pair_coeff.set(p,
                                 p + '_binder',
                                 epsilon=settings.epsilon,
                                 sigma=1.0,
                                 r_cut=False)
    # No attraction for grey bead and binders
    if (p == 'A'):
        continue
    lj_att.pair_coeff.set(p,
                          p + '_binder',
                          epsilon=settings.epsilon,
                          sigma=1.0,
                          r_cut=settings.lj_rcut_att)
    lj_rep.pair_coeff.set(p,
                          p + '_binder',
                          epsilon=settings.epsilon,
                          sigma=1.0,
                          r_cut=False)

hoomd.run_upto(settings.final_timesteps)
restart_gsd.write_restart()
