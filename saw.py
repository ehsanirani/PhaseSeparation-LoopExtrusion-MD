import numpy as np
import random


def pbc(r, D):
    Dhalf = 0.5 * D
    res = r
    system_dims = r.shape[0]
    for e in range(0, system_dims):
        if res[e] < -Dhalf[e]:
            res[e] = res[e] + D[e]
        elif res[e] >= Dhalf[e]:
            res[e] = res[e] - D[e]
    return res


def get_init_pol_xyz(N=100,
                     sigma=1.0,
                     step_mag=0.97,
                     no_backfolding_dist=1.02,
                     D=0,
                     pbc_flag=np.ones(3)):
    random.seed()
    R = 0.5 * sigma  # beads radius
    walk = np.zeros(shape=(N, 3))  # the walk chain!
    img = np.zeros(shape=(N, 3))
    # generate random chain:
    n = 0
    nmolecule = 1
    # random starting point for the chain in the box
    r1 = np.array([0.0, 0.0, 0.0])
    r2 = np.array([random.random(), random.random(), random.random()]) * D
    # store 1st monomer of chain
    # ist monomer is always in original box (image = 0)
    r2 = pbc(r2, D)  # need to be coded
    walk[n] = r2
    img[n] = np.array([0, 0, 0])
    n = n + 1
    # generate rest of monomers in the chain
    for imonomer in range(1, N):
        r0 = r1
        r1 = r2
        # random points inside sphere of unit radius
        accept = False
        while (accept == False):
            r = 2.0
            while (r > 1.0):
                rinner = np.array([
                    random.uniform(-1, 1),
                    random.uniform(-1, 1),
                    random.uniform(-1, 1)
                ])
                r = np.linalg.norm(rinner)

            accept = True
            # project point to surface of sphere of unit radius
            rsurf = rinner / r
            # create new point by scaling unit offsets by bondlength (sigma = 1.0)
            r2 = r1 + rsurf * step_mag
            # check that new point meets restriction on backfolding
            dr = r2 - r0
            r = np.linalg.norm(dr)
            if (imonomer > 2) and (r <= no_backfolding_dist):
                accept = False
            # check pbc:
            crossed_boundaries = np.abs(r2) > D * 0.5
            accept = not ((pbc_flag.astype(int) -
                           crossed_boundaries.astype(int)) == -1).any()

        # store new point
        # if delta to revious bead is large,  then increment/decrement image flag
        r2 = pbc(r2, D)  # need to be coded
        walk[n] = r2
        dr = walk[n] - walk[n - 1]
        for e in range(0, 3):
            if abs(dr[e]) < 2.0 * step_mag:
                img[n, e] = img[n - 1, e]
            elif dr[e] < 0.0:
                img[n, e] = img[n - 1, e] + 1
            elif dr[e] > 0.0:
                img[n, e] = img[n - 1, e] - 1
        n = n + 1
    return walk, img
