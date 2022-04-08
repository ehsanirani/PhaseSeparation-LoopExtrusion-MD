import saw
import numpy as np

class sbs:
    def __init__(self, hoomd, mode):
        self.hoomd = hoomd
        self.features = {}
        self.snap = None

    def init_world(self, mode='gpu', args=''):
        hoomd = self.hoomd
        hoomd.context.initialize('--mode={} {}'.format(mode, args))

    def set_restart(self, flag, restart_file):
        self.restartable = flag
        self.restart_file = restart_file

    def add_random_polymer(self, snap, L, **kargs):
        size = kargs['size']
        if size < 1:
            return -1
        old_N = snap.particles.N
        snap.particles.resize(old_N + size)
        pos, img = saw.get_init_pol_xyz(N=size, D=L, pbc_flag=kargs['pbc'])
        if 'rcm_pos' in kargs.keys():
            unwrapped_position = pos + img * L
            rcm = np.mean(unwrapped_position, axis=0)
            central_pos = unwrapped_position - rcm
            pos = central_pos
        snap.particles.position[-size:] = pos
        snap.particles.image[-size:] = img
        snap.particles.typeid[-size:] = kargs['particles_typeid']

        snap.bonds.resize(snap.bonds.N + size - 1)
        snap.bonds.typeid[-size + 1:] = kargs['bonds_typeid']
        snap.bonds.group[-size + 1:] = [[i, i + 1]
                                        for i in range(old_N, old_N + size - 1)
                                        ]

        if ('angles_typeid' in kargs.keys()):
            if kargs['angles_typeid'] >= 0:
                snap.angles.resize(snap.angles.N + size - 2)
                snap.angles.typeid[-size + 2:] = kargs['angles_typeid']
                snap.angles.group[-size + 2:] = [
                    [i, i + 1, i + 2] for i in range(old_N, old_N + size - 2)
                ]

    def add_random_particles(self, snap, L, **kargs):
        size = kargs['size']
        if size < 1: return -1
        old_N = snap.particles.N
        snap.particles.resize(old_N + size)
        snap.particles.position[-size:] = np.random.rand(size, 3) * L - 0.5 * L
        snap.particles.typeid[-size:] = kargs['particles_typeid']

    def remove_overlaps(self, group, limit=0.1, steps=5000, **fire_args):
        default_args = {'Etol': 1e-7, 'ftol': 1e-5, 'dt': 0.001}
        args = {**default_args, **fire_args}
        print(args)
        hoomd = self.hoomd
        md = self.hoomd.md
        nve = md.integrate.nve(group=group, limit=limit)
        fire = md.integrate.mode_minimize_fire(**args)
        counter = 0
        while not (fire.has_converged()):
            hoomd.run(steps)
            counter += 1
            if counter > 100: break
        nve.disable()

    def thermalize(self):  # code me!
        return 0

    def set_pairs(self, feature):  # code me!
        return 0
