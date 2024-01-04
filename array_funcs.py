import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
from scipy import signal as sp
from circ_rps import *
import warnings
warnings.filterwarnings("error")
plt.rc('font',family='serif')

def array_factor(xs, ys,k, f,theta,phi,t0=0,p0=0):
    if len(xs)!=len(ys):
        print('X and Y position array lengths do not match')
        raise ValueError
        return
    t0 = np.deg2rad(t0)
    p0 = np.deg2rad(p0)
    ArrF = np.zeros([len(theta),len(phi)],dtype=np.complex64)
    Tot = np.zeros([len(theta),len(phi)],dtype=np.complex64)
    for i in range(len(xs)):
        pld = (xs[i]*np.sin(t0)*np.cos(p0) + ys[i]*np.sin(t0)*np.sin(p0))
        steering_vector = np.exp(-1j*k*pld)
        u = xs[i]*np.outer(np.sin(theta),np.cos(phi))
        v = ys[i]*np.outer(np.sin(theta),np.sin(phi))
        this_element = np.exp(1j*k*(u+v))*steering_vector
        ArrF += this_element
        if np.shape(ArrF) != np.shape(f[1]):
            print(f'Element pattern data (shape: {np.shape(f)}) is different shape from Theta x Phi (shape: {np.shape(ArrF)}')
        else:
            Tot += this_element*f[i]

    ArrF = np.real(ArrF)
    ArrF = ArrF*(np.max(ArrF)/np.mean(ArrF))
    Tot = np.real(Tot)
    Tot = np.real(Tot)
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
                if (Arrf[idx,ph]-Arrf[90,90]) > psll:
                    psll = Arrf[idx,ph]-Arrf[90,90]
    asll = asll/N
    print(f'PSLL: {psll}')
    print(f'ASLL: {asll}')
    #cost = np.abs(des_bw - meas_bw) + psll
    cost = psll + asll
    
    return cost,peaks

def makeUVPlot(theta,phi,G,title,filename):
    u = np.outer(np.sin(theta),np.cos(phi))
    v = np.outer(np.sin(theta),np.sin(phi))
    G[u**2+v**2>1] = np.nan 

    fig,ax = plt.subplots(layout='constrained')
    CS = ax.pcolormesh(u,v,G,shading='gouraud',cmap='plasma',vmax=np.nanmax(G),vmin=-10)
    ax.set_xlabel('u',weight='bold',fontsize=14)
    ax.set_ylabel('v',weight='bold',fontsize=14)
    fig.colorbar(CS,label='Gain (dBi)')
    ax.set_title(title,weight='bold',fontsize=14)
    fig.savefig(filename, dpi=600)
    return

def Beamwidth(theta,Arrf):
    del_theta = (theta[1] - theta[0])*180/np.pi
    try:
        peak = sp.peak_widths(Arrf[:,90],[90],rel_height=0.1)[0]
    except RuntimeWarning:
        peak = 180
    return peak*del_theta

def makeArrayPlot(xs,ys,d,spirads,idstring):
    fig,ax = plt.subplots()
    xmax = 1.1*np.max(xs)+np.max(spirads)
    ymax = 1.1*np.max(ys)+np.max(spirads)
    ax.axis([-xmax, xmax, -ymax, ymax])
    circs = []
    for i in range(len(xs)):
        circs.append(Circle((xs[i],ys[i]),spirads[i],facecolor='blue',edgecolor='k'))
        ax.add_patch(circs[i])
    ax.grid()
    ax.set_aspect('equal','box')
    fig.savefig(f'results/plots/array_plot_{idstring}.png')
    plt.close()
    

if __name__ == '__main__':
    import matplotlib.pyplot as plt
    d= rps(6e9,1,3)
    xs, ys = circ_positions(d)
    xs = np.insert(xs,0,0)
    ys = np.insert(ys,0,0)
    rs = np.array([5,-10,3,20])
    #xs,ys = rotate(xs,ys,rs,d)
    spirad = 10
    spirads = np.ones(xs.size)*spirad/1000
    freq = 6e9
    lamb = 3e8/freq
    k = 2*np.pi/lamb
    theta = np.linspace(0, np.pi, 181)
    phi = np.linspace(0, 2*np.pi, 361)
    f = np.ones((xs.size,181,361))

    # broadside pattern
    ArrF,Tot = array_factor(xs,ys,k,f,theta,phi,t0=0,p0=0)
    G = 10*np.log10(ArrF)
    title = f'No Coupling, f={freq/1e9}GHz, Broadside Scan'
    filename = f'nc_f{freq/1e9}ghz_ph0th0.pdf'
    makeUVPlot(theta,phi,G,title,filename)
    plt.show()
    
    # scanned pattern
    ArrF,Tot = array_factor(xs,ys,k,f,theta,phi,t0=60,p0=30)
    G = 10*np.log10(ArrF)
    title = f'No Coupling, f={freq/1e9}GHz, Scanned to 60 degrees'
    filename = f'nc_f{freq/1e9}ghz_ph30th60.pdf'
    makeUVPlot(theta,phi,G,title,filename)

    plt.show()
    makeArrayPlot(xs,ys,d,spirads,'uniform')
    print('debug')

