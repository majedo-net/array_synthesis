import numpy as np
import matplotlib.pyplot as plt

from scipy import signal as sp
from circ_rps import *
import warnings
warnings.filterwarnings("error")

def array_factor(xs, ys,k, f,theta,phi):
    if len(xs)!=len(ys):
        print('X and Y position array lengths do not match')
        raise ValueError
        return
    ArrF = np.zeros([len(theta),len(phi)],dtype=np.complex64)
    Tot = np.zeros([len(theta),len(phi)],dtype=np.complex64)
    for i in range(len(xs)):
        r = np.sqrt(xs[i]**2 + ys[i]**2)
        xvec = r*np.outer(np.sin(theta),np.cos(phi))
        yvec = r*np.outer(np.sin(theta),np.sin(phi))
        ArrF += np.exp(1j*k*(xvec + yvec))
        if np.shape(ArrF) != np.shape(f[1]):
            print(f'Element pattern data (shape: {np.shape(f)}) is different shape from Theta x Phi (shape: {np.shape(ArrF)}')
        else:
            Tot = ArrF*f[i]
    return np.abs(ArrF),np.abs(Tot)

def BeamCost(des_bw,meas_bw,theta,phi,Arrf):
    des_bw = des_bw
    meas_bw = meas_bw
    peaks=[]
    psll = 0
    for ph in [0,15,30,45,60,75,90,105,120,135,150,165,180]:
        peaks.extend(sp.find_peaks(Arrf[:,ph],prominence=20,width=5)[0].tolist())
        for p in peaks:
            if np.abs(np.rad2deg(theta[p])) > (meas_bw/2):
                if Arrf[p,ph] > psll:
                    psll = Arrf[p,ph]
            else:
                peaks.remove(p)
    print(f'PSLL: {psll}')
    cost = np.abs(des_bw - meas_bw) + psll
    
    return cost,peaks

def Beamwidth(theta,Arrf):
    del_theta = (theta[1] - theta[0])*180/np.pi
    try:
        peak = sp.peak_widths(Arrf[:,90],[90],rel_height=0.5)[0]
    except RuntimeWarning:
        peak = 180
    return peak*del_theta

def makePatternPlots(theta,phi,AF,Tot,cost,freq,peaks=None,element=None,save=False):
    fig,[ax1,ax2,ax3] = plt.subplots(3,1)
    for ph in [0,15,30,45,60,75,90,105,120,135,150,165,180]:
        ph = int(ph)
        ax1.plot(np.rad2deg(phi),Tot[ph,:],label=f'Theta={ph}')
        ax2.plot(np.rad2deg(theta),Tot[:,ph],label=f'Phi={ph}')
        if peaks:
            ax2.plot(np.rad2deg(theta[peaks]),Tot[peaks,ph],'x')
        ax3.plot(np.rad2deg(theta),element[:,ph],label=f'phi = {ph}')
    ax1.grid(True,which='both')
    ax2.grid(True,which='both')
    ax1.legend()
    ax2.legend()
    ax1.set_title('Total')
    fig.set_size_inches(10,8)
    if save:
        fig.savefig(f'/results/freq_{freq/1e6}_cost_{int(cost)}.png')
    fig,[ax1,ax2] = plt.subplots(2,1)
    for ph in [0,15,30,45,60,75,90,105,120,135,150,165,180]:
        ph = int(ph)
        ax1.plot(np.rad2deg(phi),AF[ph,:],label=f'Theta={ph}')
        ax2.plot(np.rad2deg(theta),AF[:,ph],label=f'Phi={ph}')
    ax1.grid(True,which='both')
    ax2.grid(True,which='both')
    ax1.legend()
    ax2.legend()
    ax1.set_title('Array Factor')


if __name__ == '__main__':
    import matplotlib.pyplot as plt
    d= rps(3e9,1.5,3)
    xs, ys = circ_positions(d)
    xs = np.insert(xs,0,0)
    ys = np.insert(ys,0,0)
    rs = np.array([5,-10,3,20])
    xs,ys = rotate(xs,ys,rs,d)
    f = np.ones(len(xs),dtype=object)
    spirad = 70 
    freq = 1870e6
    lamb = 3e8/freq
    for i in np.arange(len(f)):
        pat = np.genfromtxt(f'./results/ff/farfieldspiral_rad_{spirad}_freq_{int(freq/1e6)}.csv',delimiter=',')
        pat = pat / np.amax(pat)
        f[i] = pat
    k = 2*np.pi/lamb
    theta = np.linspace(-np.pi/2, np.pi/2, 181)
    phi = np.linspace(0, np.pi, 181)
    ArrF,Tot = array_factor(xs,ys,k,f,theta,phi)

    PHI, TH= np.meshgrid(phi,theta)
    des_bw = 40 # desired beamwidth +/- degrees
    meas_bw = Beamwidth(theta,ArrF)
    
    cost,peaks=BeamCost(des_bw,meas_bw,theta,phi,Tot)
    fig,ax = plt.subplots()
    ax.plot(xs,ys,'o')
    ax.grid(True,'both')
    
    print(f'freq: {freq} \t measBW: {meas_bw} \t cost: {cost}')
    makePatternPlots(theta,phi,ArrF,Tot,peaks,pat)    

    plt.show()
