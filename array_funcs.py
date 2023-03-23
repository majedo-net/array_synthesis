import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
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
        u = r*np.outer(np.sin(theta),np.cos(phi))
        v = r*np.outer(np.sin(theta),np.sin(phi))
        ArrF += np.exp(1j*k*(u + v))
        if np.shape(ArrF) != np.shape(f[1]):
            print(f'Element pattern data (shape: {np.shape(f)}) is different shape from Theta x Phi (shape: {np.shape(ArrF)}')
        else:
            Tot = ArrF*f[i]
    return np.abs(ArrF),np.abs(Tot)

def BeamCost(des_bw,meas_bw,theta,phi,Arrf):
    des_bw = des_bw
    meas_bw = meas_bw
    print(meas_bw)
    peaks=[]
    psll = 0 # peak sidelobe level
    asll = 0 # average sidelobe level
    N = 1
    for ph in [0,15,30,45,60,75,90,105,120,135,150,165,180]:
        for idx,th in enumerate(theta):
            if np.abs(np.rad2deg(th)) > meas_bw/2:
                asll += Arrf[idx,ph]
                N= N+1
                if Arrf[idx,ph] > psll:
                    psll = Arrf[idx,ph]-Arrf[90,90]
    asll = asll/N
    print(f'PSLL: {psll}')
    print(f'ASLL: {asll}')
    #cost = np.abs(des_bw - meas_bw) + psll
    cost = psll + asll
    
    return cost,peaks

def Beamwidth(theta,Arrf):
    del_theta = (theta[1] - theta[0])*180/np.pi
    try:
        peak = sp.peak_widths(Arrf[:,90],[90],rel_height=0.5)[0]
    except RuntimeWarning:
        peak = 180
    return peak*del_theta

def makeArrayPlot(xs,ys,d,spirads,idstring):
    fig,ax = plt.subplots()
    xmax = 1.1*np.max(xs)+np.max(spirads)
    ymax = 1.1*np.max(ys)+np.max(spirads)
    ax.axis([-xmax, xmax, -ymax, ymax])
    dind = 1
    circs = []
    for i in range(len(xs)):
        circs.append(Circle((xs[i],ys[i]),spirads[i],facecolor='blue',edgecolor='k'))
        ax.add_patch(circs[i])
    ax.grid()
    ax.set_aspect('equal','box')
    fig.savefig(f'/results/plots/array_plot_{idstring}.png')
    plt.close()
    

def makePatternPlots(theta,phi,AF,Tot,cost,freq,idstring,peaks=None,element=None,save=False):
    fig,[ax1,ax2] = plt.subplots(2,1)
    for ph in [0,15,30,45,60,75,90,105,120,135,150,165,180]:
        ph = int(ph)
        ax1.plot(np.rad2deg(phi),Tot[ph,:],label=f'Theta={ph}')
        ax2.plot(np.rad2deg(theta),Tot[:,ph],label=f'Phi={ph}')
        if peaks:
            ax2.plot(np.rad2deg(theta[peaks]),Tot[peaks,ph],'x')
    ax1.grid(True,which='both')
    ax2.grid(True,which='both')
    ax1.legend()
    ax2.legend()
    ax1.set_title('Total')
    fig.set_size_inches(10,8)
    if save:
        fig.savefig(f'/results/plots/id_{idstring}_freq_{int(freq/1e6)}_cost_{int(cost)}.png')

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
