# -*- coding: utf-8 -*-
"""
Created on Tue Oct 18 15:19:47 2022

@author: Matt Dodd
"""
import numpy as np
from array_funcs import *

def circ_count(len_d, count=0):
    count = 6+(2*len_d) + count
    if len_d <=0:
        return count
    len_d = len_d-1
    return circ_count(len_d,count)

def rps(fmax,r,N):
    d0 = (3e8/fmax)/2
    if r >=1:
        xi = 1
    else:
        xi = 1/((N**r) - (N-1)**r)
    d= np.zeros(N+1)
    for n in range(N+1):
        d[n] = d0*xi*(n+1)**r 
    
    return d

def rotate(xs,ys,theta,d):
    if len(theta) > len(d):
        print('Received more rotation angles than layers...')
        return
    theta = theta * np.pi / 180
    rotmat = np.array([[np.cos(theta), -1*np.sin(theta)],
                       [np.sin(theta), np.cos(theta)]])
    for n,t in enumerate(theta):
        if n==0: b=0
        else: b = circ_count(n-1)
        if n>0: b+=1
        e = circ_count(n) +1
        [xs[b:e],ys[b:e]] = np.matmul(rotmat[:,:,n],[xs[b:e],ys[b:e]])
    return xs,ys
                       
    

def circ_positions(d):
    xs = [] 
    ys = [] 
    for n,dn in enumerate(d):
        n_layer = (2*n) +6
        wd = (2*np.pi/n_layer)
        for w in np.arange(0,(2*np.pi),wd):
            xs.append(dn*np.cos(w))
            ys.append(dn*np.sin(w))

    xs = np.insert(xs,0,0)
    ys = np.insert(ys,0,0)
    return xs,ys
if __name__ == '__main__':
    import matplotlib.pyplot as plt
    d= rps(3e9,0.8,3)
    xs, ys = circ_positions(d)
    rs = np.array([20,-10.3,-1.58,12.58])
    #xs,ys = rotate(xs,ys,rs,d)
    plt.scatter(xs,ys)
    plt.grid()
    ax = plt.gca()
    ax.set_aspect('equal','box')
    plt.show()
