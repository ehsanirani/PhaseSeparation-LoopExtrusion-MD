import numpy as np


class LEagents:
    """ write some docstring """
    def __init__(self,
                 system=None,
                 stop_probs=None,
                 N_le=0,
                 p_off=0.0125,
                 update_period=1,
                 move_dist_thr=-1,
                 le_typeid=1,
                 le_type='loop_extrusion',
                 passing_p=1.0,
                 CTCF_type='all'):
        self.system = system
        self.box = system.box
        self.L = np.array([self.box.Lx, self.box.Ly, self.box.Lz])
        self.Lhalf = 0.5 * self.L
        if (CTCF_type == 'all'):
            self.stop_probs = stop_probs
        elif (CTCF_type == 'random'):
            self.stop_probs = (np.random.rand(*stop_probs.shape) <
                               stop_probs).astype(float) * passing_p
        self.p_off = p_off
        self.update_period = update_period
        self.move_dist_thr = move_dist_thr
        self.le_typeid = le_typeid
        self.le_type = le_type
        self.N_le = N_le
        self.new_le_N = 0
        self.remove_le_N = 0
        self.passing_p = passing_p

    def replace_loop_extrusion_factors(self, snap):
        """
        Replace current LE factors with new ones, according to given rates.
        """
        # Removing part
        system = self.system
        # remove_le_N = self.remove_le_N # old way. Now it will depend only on p_off
        s = np.random.sample(self.N_le) < self.p_off
        remove_le_N = s.sum()
        self.remove_le_N = remove_le_N
        if remove_le_N == 0:
            return  # No LE is selected to dissociate
        tags = []
        for b in self.system.bonds:
            if b.typeid != self.le_typeid:
                continue
            tags.append(b.tag)
        if len(tags) > 0:
            to_be_deleted = np.random.choice(tags, remove_le_N, replace=False)
            if len(to_be_deleted) > 0:
                for t0 in to_be_deleted:
                    system.bonds.remove(t0)
            self.N_le = len(tags) - len(to_be_deleted)

        # Creating part
        system = self.system
        N_beads = len(self.stop_probs)
        beads_list = np.arange(N_beads)
        new_le_N = remove_le_N
        self.new_le_N = new_le_N
        tags = []
        for b in self.system.bonds:
            if b.typeid != self.le_typeid:
                continue
            tags.append(b.tag)
        ab_list = np.zeros((len(tags), 2))
        for i, t in enumerate(tags):
            b = system.bonds.get(t)
            ab_list[i, 0] = b.a
            ab_list[i, 1] = b.b
        # p = np.ones(N_beads)
        options = beads_list[:-1]
        if len(ab_list) > 0:
            le_lengths = np.diff(ab_list)[:, 0]
            not_permitted_a = ab_list[le_lengths == 1][:, 0].tolist()
            if len(not_permitted_a) > 0:
                options = [
                    x for x in beads_list[:-1] if x not in not_permitted_a
                ]
        new_a = np.random.choice(options, new_le_N, replace=False)
        for a0 in new_a:
            self.system.bonds.add(self.le_type, a0, a0 + 1)

        self.N_le = len(tags) + len(new_a)

    def update_loop_extrusion2(self, timestep):
        def get_le_neighbs(x, direction, up_lim):
            res = x + direction
            res[res < 0] = 0
            res[res >= up_lim] = up_lim - 1
            return res

        system = self.system
        sl_typeid = self.le_typeid
        stop_probs = self.stop_probs.copy()

        snap = system.take_snapshot(all=True)
        loop_extrusions_ids = snap.bonds.typeid == sl_typeid
        loop_extrusions = snap.bonds.group[loop_extrusions_ids].astype(int)
        N_loop_extrusions = loop_extrusions.shape[0]
        if loop_extrusions.shape[0] == 0:
            return

        group1 = loop_extrusions.copy()

        # To avoid sli passing each others
        neighbs = get_le_neighbs(loop_extrusions[:, 0], -1,
                                 stop_probs.shape[0])
        stop_probs[neighbs, 1] = np.max([
            stop_probs[neighbs, 1],
            np.ones(len(loop_extrusions)) * self.passing_p
        ],
                                        axis=0)
        neighbs = get_le_neighbs(loop_extrusions[:, 1], +1,
                                 stop_probs.shape[0])
        stop_probs[neighbs, 0] = np.max([
            stop_probs[neighbs, 0],
            np.ones(len(loop_extrusions)) * self.passing_p
        ],
                                        axis=0)

        neighbs = get_le_neighbs(loop_extrusions[:, 0], +1,
                                 stop_probs.shape[0])
        stop_probs[neighbs, 0] = np.max([
            stop_probs[neighbs, 0],
            np.ones(len(loop_extrusions)) * self.passing_p
        ],
                                        axis=0)
        neighbs = get_le_neighbs(loop_extrusions[:, 1], -1,
                                 stop_probs.shape[0])
        stop_probs[neighbs, 1] = np.max([
            stop_probs[neighbs, 1],
            np.ones(len(loop_extrusions)) * self.passing_p
        ],
                                        axis=0)

        ######

        random_probs = np.random.rand(loop_extrusions.shape[0], 2)

        accept_moves = np.zeros_like(random_probs).astype(bool)
        accept_moves[:, 0] = stop_probs[loop_extrusions[:, 0],
                                        0] < random_probs[:, 0]
        accept_moves[:, 1] = stop_probs[loop_extrusions[:, 1],
                                        1] < random_probs[:, 1]

        check_dist_ids = np.logical_or(accept_moves[:, 0], accept_moves[:, 1])
        dist_small = np.zeros_like(check_dist_ids).astype(bool)

        group1[check_dist_ids, 0] -= accept_moves[check_dist_ids, 0]  # go down
        group1[check_dist_ids, 1] += accept_moves[check_dist_ids, 1]  # fo up
        #group1[group1>=snap.particles.N] = snap.particles.N-1
        #group1[group1<0] = 0

        # check if we have a distance threshold to extrude
        # dist. thresh. -1 means no threshold!
        if self.move_dist_thr != -1:
            dr0 = self.dist_arr(snap, group1[check_dist_ids])
            dist_small[check_dist_ids] = dr0 < self.move_dist_thr
        else:
            dist_small[:] = True

        not_update_ids = np.logical_not(dist_small)

        # reverse go down
        group1[not_update_ids, 0] += accept_moves[not_update_ids, 0]
        # reverse fo up
        group1[not_update_ids, 1] -= accept_moves[not_update_ids, 1]

        snap.bonds.group[loop_extrusions_ids] = group1
        system.restore_snapshot(snap)
        ##
        # replace LEagents
        self.replace_loop_extrusion_factors(snap)
        ######
