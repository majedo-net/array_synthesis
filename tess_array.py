import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
from scipy import signal as sp
from scipy.spatial import KDTree
from circ_rps import *
import warnings,os,tempfile
warnings.filterwarnings("error")
plt.rc('font',family='serif')

import spiral
'''
from CSXCAD  import ContinuousStructure
from openEMS import openEMS
from openEMS.physical_constants import *
'''

def getNearestNeighbors(xs,ys,idx,neighbors):
    points = list(zip(xs.ravel(),ys.ravel()))
    tree = KDTree(points)
    ds,idxes = tree.query([xs[idx],ys[idx]],k=(neighbors+1))
    txs = xs[idxes[1:]]
    tys = ys[idxes[1:]]
    return txs,tys


if __name__ == '__main__':
    import matplotlib.pyplot as plt
    d= rps(6e9,1,3)
    xs, ys = circ_positions(d)
    rs = np.array([5,-10,3,20])
    #xs,ys = rotate(xs,ys,rs,d)
    spirad = 10
    spirads = np.ones(xs.size)*spirad
    freq = 5e9
    lamb = 3e8/freq
    k = 2*np.pi/lamb
    theta = np.linspace(0, np.pi, 181)
    phi = np.linspace(0, 2*np.pi, 361)
    fs = np.zeros([xs.size, theta.size, phi.size])
    getNearestNeighbors(xs,ys,0,8)