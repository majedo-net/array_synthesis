import numpy as np
import matplotlib.pyplot as plt
import multiprocessing as mp
from scipy.spatial import KDTree
from circ_rps import *
import warnings,datetime
#warnings.filterwarnings("error")
plt.rc('font',family='serif')

import eqspiral as eqsp
import patch as psim
from array_funcs import *

def getNearestNeighbors(xs,ys,idx,neighbors):
    points = list(zip(xs.ravel(),ys.ravel()))
    tree = KDTree(points)
    ds,idxes = tree.query([xs[idx],ys[idx]],k=(neighbors+1))
    txs = xs[idxes]
    tys = ys[idxes]
    return txs,tys,idxes

def getElementParams(xs,ys,eid,freq,hs,L,W,theta,phi):
    # prepare element parameters for a multi processing worker
    txs,tys,neighbor_indexes = getNearestNeighbors(xs,ys,eid,6)
    centers = np.vstack((txs,tys)).T*1000
    centers = centers - centers[0,:]
    return [freq,hs,centers,neighbor_indexes,L,W,theta,phi,eid]

def parallelEmbeddedWorker(element):
    freq,hs,centers,neighbor_indexes,L,W,theta,phi,eid = element

    En,s11f,s11_db,sfreq,zin,sn1 = psim.SimulateEmbeddedFarfield(freq,hs,centers,L,W,theta,phi,eid=eid)
    print(len(neighbor_indexes))
    print(type(neighbor_indexes))
    neighbor_indexes = np.asarray(neighbor_indexes,dtype=np.int32)
    np.savetxt(f'/results/ff_{eid}.txt',En)
    np.savetxt(f'/results/s11_{eid}.txt',(sfreq,s11_db,np.real(zin),np.imag(zin)))
    np.savetxt(f'/results/sn{eid}.txt',(neighbor_indexes,sn1[:,151]))

    return


'''
Single threaded version
=======================
'''
if __name__ == '__main__':
    import matplotlib.pyplot as plt
    with open('/results/times.txt',mode='a') as f:
        f.write(f'start: {datetime.datetime.now()}\n')
    #d= rps(6e9,1,3)
        d= rps(6e9,1,1)
    xs, ys = circ_positions(d)
    rs = np.array([5,-10,3,20])
    #xs,ys = rotate(xs,ys,rs,d)
    L = 10.2
    W = 15.5
    hs =1.5 
    freqs = [6e9]
    for fdx in range(len(freqs)):
        freq = freqs[fdx]
        lamb = 3e8/freq
        k = 2*np.pi/lamb
        theta = np.linspace(0, np.pi, 181)
        phi = np.linspace(0, 2*np.pi, 361)
        fs = np.zeros([xs.size, theta.size, phi.size])
        # -----------------------------#
        # Get embedded Element Patterns
        # -----------------------------#
        print('=============================================')
        print('Starting Embedded Element Pattern Simulation')
        print('=============================================')

        centers = np.vstack((xs,ys)).T*1000
        smn = psim.SimulateFullArray(freq,hs,centers,L,W,theta,phi)
        np.savetxt('/results/fullarray_smn.txt',smn)
    with open('/results/times.txt',mode='a') as f:
        f.write(f'end: {datetime.datetime.now()}')

