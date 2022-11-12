import numpy as np

from scipy import signal as sp
from hex_rps import *

def array_factor(xs, ys,k, f,theta,phi):
    if len(xs)!=len(ys):
        print('X and Y position array lengths do not match')
        raise ValueError
        return
    ArrF = np.zeros([len(theta),len(phi)],dtype=np.complex64)
    for i in range(len(xs)):
        xvec = xs[i]*np.outer(np.cos(theta),np.cos(phi))
        yvec = ys[i]*np.outer(np.sin(theta),np.sin(phi))
        ArrF += np.exp(1j*k*(xvec + yvec))
    
    return np.abs(ArrF)

def BeamCost(des_bw,meas_bw,sll,theta,phi,Arrf):
    des_bw = des_bw*np.pi/180
    meas_bw = meas_bw*np.pi/180
    slIdxsth = np.argwhere(np.abs(theta) > (meas_bw)/2).T
    slIdxsph = np.argwhere(np.abs(phi-np.pi/2) < (meas_bw)/2)
    mlIdxsth = np.argwhere(np.abs(theta) < (meas_bw)/2).T
    mlIdxsph = np.argwhere(np.abs(phi-np.pi/2) < (meas_bw)/2)
    sll = 10**(-1*sll/10)
    cost = np.max(Arrf[slIdxsth,slIdxsph])*sll - np.max(Arrf[mlIdxsth,mlIdxsph]) + np.abs(des_bw - meas_bw)*180/np.pi
    
    return cost

def Beamwidth(theta,Arrf):
    del_theta = (theta[1] - theta[0])*180/np.pi
    peak1 = sp.peak_widths(Arrf[90,:],[90],rel_height=0.2)[0]
    peak2 = sp.peak_widths(Arrf[:,90],[90],rel_height=0.2)[0]
    return np.max([peak1,peak2])*del_theta


if __name__ == '__main__':
    import matplotlib.pyplot as plt
    d= rps(3e9,0.89,3)
    xs, ys = hex_positions(d)
    xs = np.insert(xs,0,0)
    ys = np.insert(ys,0,0)
    rs = np.array([1.2,-10.3,-1.58,12.58])
    xs,ys = rotate(xs,ys,rs,d)
    f = np.ones(len(xs))
    freq = 400e6
    lamb = 3e8/freq
    k = 2*np.pi/lamb
    theta = np.linspace(-np.pi/2, np.pi/2, 181)
    phi = np.linspace(0, np.pi, 181)
    ArrF = array_factor(xs,ys,k,f,theta,phi)

    PHI, TH= np.meshgrid(phi,theta)
    des_bw = 40 # desired beamwidth +/- degrees
    sll = -5 # desired sidelobe level (used for weighting)
    meas_bw = Beamwidth(theta,ArrF)
    
    cost=BeamCost(des_bw,meas_bw,sll,theta,phi,ArrF)
    
    print(f'freq: {freq} \t measBW: {meas_bw} \t cost: {cost}')
    X = ArrF * np.cos(PHI) * np.cos(TH)
    Y = ArrF * np.sin(PHI) * np.sin(TH)
    Z = ArrF * np.cos(TH)
    fig = plt.figure()
    # ax = fig.add_subplot(1,1,1)
    # ax.matshow(ArrF)
    fig2 = plt.figure()
    ax2 = fig2.add_subplot(1,1,1)
    ax2.scatter(xs,ys)
    ax2.grid()
    ax = fig.add_subplot(1,1,1, projection='3d')
    plot = ax.plot_surface(
        X,Y,Z, cmap=plt.get_cmap('jet'),  
        linewidth=0, antialiased=False, alpha=0.5)

    plt.show()