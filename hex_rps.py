# -*- coding: utf-8 -*-
"""
Created on Tue Oct 18 15:19:47 2022

@author: Matt Dodd
"""
import numpy as np
from array_funcs import *

def hex_count(len_d, elements=0):
    elements = len_d*6 + elements
    if len_d ==0:
        return elements
    len_d = len_d-1
    return hex_count(len_d,elements)

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
        b = hex_count(n) + 1
        e = hex_count(n+1) + 1
        [xs[b:e],ys[b:e]] = np.matmul(rotmat[:,:,n],[xs[b:e],ys[b:e]])
    return xs,ys
                       
    

def hex_positions(d):
    xs = [None] * hex_count(len(d))
    ys = [None] * hex_count(len(d))
    for n,dn in enumerate(d):
        #dnsq = np.sqrt((dn**2)/2)
        dnsq = dn/2
        if n==0:
            xs[0] = dn 
            ys[0] = 0
            xs[1] = dnsq
            ys[1] = dnsq
            xs[2] = -dnsq
            ys[2] = dnsq
            xs[3] = -dn
            ys[3] = 0
            xs[4] = -dnsq
            ys[4] = -dnsq
            xs[5] = dnsq
            ys[5] = -dnsq
        else:
            xs[hex_count(n)] = dn 
            ys[hex_count(n)] = 0
            dn =  d[n]/(n+1)
            dnsq = d[n]/(n+1)/2
            for i in range(hex_count(n)+1,hex_count(n)+(n+2)):
                xs[i] = xs[i-1] - dnsq
                ys[i] = ys[i-1] + dnsq
            for i in range(hex_count(n)+(n+2),hex_count(n)+(2*n+3)):
                xs[i] = xs[i-1] - dn
                ys[i] = ys[i-1] 
            for i in range(hex_count(n)+(2*n+3),hex_count(n)+(3*n+4)):
                xs[i] = xs[i-1] - dnsq
                ys[i] = ys[i-1] - dnsq
            for i in range(hex_count(n)+(3*n+4),hex_count(n)+(4*n+5)):
                xs[i] = xs[i-1] + dnsq
                ys[i] = ys[i-1] - dnsq
            for i in range(hex_count(n)+(4*n+5),hex_count(n)+(5*n+6)):
                xs[i] = xs[i-1] + dn
                ys[i] = ys[i-1]
            for i in range(hex_count(n)+(5*n+6),hex_count(n+1)):
                xs[i] = xs[i-1] + dnsq
                ys[i] = ys[i-1] + dnsq
    
    xs = np.insert(xs,0,0)
    ys = np.insert(ys,0,0)
    return xs,ys
if __name__ == '__main__':
    import matplotlib.pyplot as plt
    d= rps(3e9,1.8,3)
    xs, ys = hex_positions(d)
    rs = np.array([2,-10.3,-1.58,12.58])
    #xs,ys = rotate(xs,ys,rs,d)
    plt.scatter(xs,ys)
    plt.grid()
    ax = plt.gca()
    ax.set_aspect('equal','box')
    plt.show()
