import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
from scipy import signal as sp
from scipy.spatial import KDTree
from circ_rps import *
import warnings,os,tempfile
warnings.filterwarnings("error")
plt.rc('font',family='serif')

import spiral.spiral as eqsp
from array_funcs import *

def getNearestNeighbors(xs,ys,idx,neighbors):
    points = list(zip(xs.ravel(),ys.ravel()))
    tree = KDTree(points)
    ds,idxes = tree.query([xs[idx],ys[idx]],k=(neighbors+1))
    txs = xs[idxes]
    tys = ys[idxes]
    return txs,tys


if __name__ == '__main__':
    import matplotlib.pyplot as plt
    d= rps(6e9,1,3)
    xs, ys = circ_positions(d)
    rs = np.array([5,-10,3,20])
    #xs,ys = rotate(xs,ys,rs,d)
    spirad = 10
    spirads = np.ones(xs.size)*spirad
    hs = 0.1
    h = 5
    freq = 5e9
    lamb = 3e8/freq
    k = 2*np.pi/lamb
    theta = np.linspace(0, np.pi, 181)
    phi = np.linspace(0, 2*np.pi, 361)
    fs = np.zeros([xs.size, theta.size, phi.size])
    # -----------------------------#
    # Get embedded Element Patterns
    # -----------------------------#
    for eid in range(xs.size):
        txs,tys = getNearestNeighbors(xs,ys,eid,8)
        centers = np.vstack((txs,tys)).T
        centers = centers - centers[0,:]
        radii = np.ones((txs.size,2))
        radii[:,0] = 0.5
        radii[:,1] = spirad
        En,s11,sfreq = eqsp.SimulateEmbeddedFarfield(freq,hs,h,centers,radii,theta,phi)
        np.savetxt(f'/results/s11_{eid}.txt',(sfreq,s11))
        fs[eid] = En

    # broadside pattern
    ArrF,Tot = array_factor(xs,ys,k,f,theta,phi,t0=0,p0=0)
    G = 10*np.log10(ArrF)
    title = f'No Coupling, f={freq/1e9}GHz, Broadside Scan'
    filename = f'/results/nc_f{freq/1e9}ghz_ph0th0.pdf'
    makeUVPlot(theta,phi,G,title,filename)
    plt.show()
    
    # scanned pattern
    ArrF,Tot = array_factor(xs,ys,k,f,theta,phi,t0=60,p0=30)
    G = 10*np.log10(ArrF)
    title = f'No Coupling, f={freq/1e9}GHz, Scanned to 60 degrees'
    filename = f'/results/nc_f{freq/1e9}ghz_ph30th60.pdf'
    makeUVPlot(theta,phi,G,title,filename)
