import numpy as np
from scipy.spatial import KDTree
from patchAntenna import PatchAntenna

class Array():
    '''
    Class to represent an antenna array

    '''

    def __init__(self,_freqs):
        self.freqs = _freqs
        self.elements = []
        return

    def circ_count(self,len_d, count=0):
        count = 6+(2*len_d) + count
        if len_d <=0:
            return count
        len_d = len_d-1
        return circ_count(len_d,count)

    def rps(self,fmax,r,N):
        d0 = (3e8/fmax)/2
        if r >=1:
            xi = 1
        else:
            xi = 1/((N**r) - (N-1)**r)
        d= np.zeros(N+1)
        for n in range(N+1):
            d[n] = d0*xi*(n+1)**r 
        return d

    def circ_positions(self,d):
        self.xs = [] 
        self.ys = [] 
        for n,dn in enumerate(d):
            n_layer = (2*n) +6
            wd = (2*np.pi/n_layer)
            for w in np.arange(0,(2*np.pi),wd):
                self.xs.append(dn*np.cos(w))
                self.ys.append(dn*np.sin(w))

        self.xs = np.insert(self.xs,0,0)
        self.ys = np.insert(self.ys,0,0)
    
    def rotate_positions(self,rot_angles):
        theta = rot_angles * np.pi / 180
        rotmat = np.array([[np.cos(theta), -1*np.sin(theta)],
                       [np.sin(theta), np.cos(theta)]])
        for n,t in enumerate(theta):
            if n==0: b=0
            else: b = self.circ_count(n-1)
            if n>0: b+=1
            e = self.circ_count(n) +1
            [self.xs[b:e],self.ys[b:e]] = np.matmul(rotmat[:,:,n],[self.xs[b:e],self.ys[b:e]])

    def generateRPSPositions(self,r=1,Nrps=3,rotations=None):
        self.circ_positions(self.rps(r,Nrps))
        self.N = len(self.xs)
        if rotations is not None:
            self.rotate_positions(rotations)
        self.xmax = np.max(np.abs(self.xs))
        self.ymax = np.max(np.abs(self.ys))

    def initPatchElements(self):
        for id in range(len(self.xs)):
            self.elements.append(PatchAntenna(
                id,
                L_=10.2,
                W_=15.5,
                x_=self.xs[id],
                y_=self.ys[id],
                hs_=1.5,
                epsr_=4.2))





