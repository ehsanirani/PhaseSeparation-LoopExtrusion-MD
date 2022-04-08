#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jun  2 08:22:54 2018

@author: ehsan
"""

import gsd, gsd.hoomd, argparse, numpy as np
import scipy.spatial.distance as spdist
import gzip
import h5py

parser = argparse.ArgumentParser()
parser.add_argument('-i', '--input-gsd', type=str, dest='input_file', help='the input gsd file')
parser.add_argument('-o', '--out', default='map', type=str, dest='fout', help='the output HDF file')
parser.add_argument('-k', '--key', default='key', type=str, dest='key', help='the key to store the results with it in the HDF file')
parser.add_argument('--particles', default=0, type=int, dest='particles', help='either 0 (all) or the number first particles which should be used. default is 0' )
parser.add_argument('--t1', type=float, default=0.0, dest='t1', help='relative start point for sampling (between 0. and 1.)')
parser.add_argument('--t2', type=float, default=1.0, dest='t2', help='relative end point for sampling (between 0. and 1.)')
parser.add_argument('--step', type=float, default=.01, dest='step', help='relative sampling step/gap size (between 0. and 1.). give 0 to use all snapshots between t1 and t2.')
parser.add_argument('--contact-thr', type=float, default=2.0, dest='contact_thr', help='distance threshold to call a contact')
parser.add_argument('--dist', type=bool, default=False, dest='dist', help='if pass True, creates the distance map as well')
parser.add_argument('--triu', type=bool, default=False, dest='triu', help='write only the upper triangle of the symetric matrices')
parser.add_argument('--gz', type=bool, default=False, dest='gz', help='if pass True, compress the output')
args = parser.parse_args()
input_file = args.input_file
#input_file = '/home/ehsan/winterfell/simulations/hoomd-saw-gen/saw-gsds/Nc5test1.gsd'
t = gsd.hoomd.open(input_file, 'rb')
#%%
def unwrapp(pos, L):
    Lhalf = 0.5*L
    res = np.zeros_like(pos)
    img = np.zeros_like(pos)
    res[0] = pos[0]
    for i in range(1,len(pos)):
        res[i] = pos[i]
        dr = res[i] - res[i-1]
        for ri in range(3):
            if dr[ri] > Lhalf:
                res[i,ri] -= L
                img[i, ri] -=1
            elif dr[ri] <-Lhalf:
                res[i,ri] += L
                img[i, ri] += 1
    return res, img
# Print iterations progress
def printProgressBar (iteration, total, prefix = '', suffix = '', decimals = 1, length = 100, fill = '█'):
    """
    Call in a loop to create terminal progress bar
    @params:
        iteration   - Required  : current iteration (Int)
        total       - Required  : total iterations (Int)
        prefix      - Optional  : prefix string (Str)
        suffix      - Optional  : suffix string (Str)
        decimals    - Optional  : positive number of decimals in percent complete (Int)
        length      - Optional  : character length of bar (Int)
        fill        - Optional  : bar fill character (Str)
    """
    percent = ("{0:." + str(decimals) + "f}").format(100 * (iteration / float(total)))
    filledLength = int(length * iteration // total)
    bar = fill * filledLength + '-' * (length - filledLength)
    print('\r%s |%s| %s%% %s' % (prefix, bar, percent, suffix), end = '\r')
    # Print New Line on Complete
    if iteration == total:
        print()

def calc_dist_cont_map_opt (pos, thr):
    N = len(pos)
    thr2 = thr*thr
    dist_map = np.zeros(shape=(N,N))
    cont_map = np.zeros(shape=(N,N))
    # generate the dr matrix
    dr_mat = np.zeros(shape=(N, N, 3))
    for i in range(N):
        dr_mat[i] = pos[i] - pos
    # calculate the squared distance matrix
    for i in range(N):
        dist_map[i] = np.einsum('ij,ji->i', dr_mat[i], dr_mat[i].T)
    dist_map = np.sqrt(dist_map)
    contacts_id = dist_map <= thr
    cont_map[contacts_id] = 1

    del dr_mat
    return dist_map, cont_map

########
# main
########
#ind = t[0].particles.charge==1
N = args.particles
if N==0:
    N = t[0].particles.N
distance_map = np.zeros(shape=(N,N))
distance0_map = np.zeros(shape=(N,N))
contact_map = np.zeros(shape=(N,N))
contact0_map = np.zeros(shape=(N,N))
L = t[0].configuration.box[0]
Lhalf = 0.5*L
l = len(t)
i=0
t1 = int(args.t1*l)
t2 = int(args.t2*l)
step = int(args.step*l)
if step < 1:
    print('step := 1')
    step=1
print('{}: {}:{}:{}'.format(l, t1,t2,step))
printProgressBar(i, t2-t1, prefix = 'Progress:', suffix = 'Complete', length = 20, fill=r'|')
for snap in t[t1:t2:step]:
    pos0 = snap.particles.position[:N]
    pos_unwrapp, img = unwrapp( pos0, L)
    #### create the output gsd file
    #R2 = calc_R2 (pos_unwrapp)
    #Rg2 = calc_Rg2 (pos_unwrapp)
    #distance_map0, contact_map0 = calc_dist_cont_map_opt(pos_unwrapp, args.contact_thr)
    ######
    distance_map0 = spdist.pdist(pos_unwrapp)
    contact_ids = distance_map0<=args.contact_thr
    contact_map0 = np.zeros_like(distance_map0)
    contact_map0[contact_ids] = 1
    contact_map0 = spdist.squareform(contact_map0)
    distance_map0 = spdist.squareform(distance_map0)
    ######
    distance_map = distance_map + distance_map0
    contact_map = contact_map + contact_map0
    printProgressBar(i*step, (t2-t1), prefix = 'Progress:', suffix = 'Complete', length = 20, fill=r'|')
    i += 1
printProgressBar(t2-t1, (t2-t1), prefix = 'Progress:', suffix = 'Complete', length = 20, fill=r'|')
distance_map = distance_map / i
contact_map = contact_map / i
#print(i)
## writing the output to hdf5 file
store = h5py.File(args.fout, 'a')

if args.dist:
    if args.triu:
        triu_ids = np.triu_indices_from(distance_map)
        temp_arr = distance_map[triu_ids]
        key = r'distance/triu/'+args.key
        store.create_dataset(key, data=temp_arr)
        #np.savetxt(args.fout+'-triu.distance', temp_arr, fmt='%.5g')
    else:
        key = r'distance/full/'+args.key
        store.create_dataset(key, data=distance_map)
        #np.savetxt(args.fout+'.distance', distance_map, fmt='%.5g')
if args.triu:
    triu_ids = np.triu_indices_from(distance_map)
    temp_arr = contact_map[triu_ids]
    key = r'contact/thr{}/triu/'.format(args.contact_thr) + args.key
    store.create_dataset(key, data=temp_arr)
    #np.savetxt(args.fout+'-triu.contact', temp_arr, fmt='%.5g')
else:
    key = r'contact/thr{}/full/'.format(args.contact_thr) + args.key
    store.create_dataset(key, data=contact_map)
    #np.savetxt(args.fout+'.contact', contact_map, fmt='%.5g')
#%%
store.close()
