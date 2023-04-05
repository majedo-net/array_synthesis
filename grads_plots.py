import numpy as np
from array_funcs import *
from circ_rps import *
import matplotlib.pyplot as plt
from matplotlib.patches import Circle

plt.rc('font', family='serif')

if __name__ == '__main__':
    # uniform array plot and cost
    d = [0.065, 0.065*2, 0.065*3,0.065*4] # lambda/2 at 2.3ghz
    spirads = 0.96*0.065/2 # outer radius of spiral
    xs,ys = circ_positions(d)
    spirads = np.ones_like(xs)*spirads
    low_pattern = np.genfromtxt('spiral/farfieldspiral_rad_31_freq_400.csv',delimiter=',')
    low_pattern = low_pattern/np.amax(low_pattern)
    high_pattern = np.genfromtxt('spiral/farfieldspiral_rad_31_freq_4200.csv',delimiter=',')
    high_pattern = high_pattern/np.amax(high_pattern)

    low_patt_array = []
    high_patt_array= []
    for i in range(len(spirads)):
        low_patt_array.append(low_pattern)
        high_patt_array.append(high_pattern)
    theta = np.linspace(-np.pi/2, np.pi/2, 181)
    phi = np.linspace(0, np.pi, 181)
    lowk = 2*np.pi*400e6/3e8
    highk = 2*np.pi*2300e6/3e8
    lowAF,lowTot = scan_array_factor(xs,ys,lowk,low_patt_array,theta,phi,np.deg2rad(60),0)
    highAF,highTot = scan_array_factor(xs,ys,highk,high_patt_array,theta,phi,np.deg2rad(60),0)
    bw = Beamwidth(theta,lowTot)
    lowcost,_=BeamCost(1,bw,theta,phi,lowTot)
    bw = Beamwidth(theta,highTot)
    highcost,_=BeamCost(1,bw,theta,phi,highTot)
    cost = lowcost+highcost
    print(f'Cost: {cost}')
## Plot starting positions of uniform array
    fig,ax = plt.subplots()
    xmax = 1.1*np.max(xs)+np.max(spirads)
    ymax = 1.1*np.max(ys)+np.max(spirads)
    ax.axis([-xmax, xmax, -ymax, ymax])
    dind = 1
    circs = []
    for i in range(len(xs)):
        circs.append(Circle((xs[i],ys[i]),spirads[i],facecolor='blue',edgecolor='k'))
        ax.add_patch(circs[i])
    ax.grid(which='both')
    ax.set_aspect('equal','box')
    ax.set_xlabel('Position (m)',fontsize=14)
    ax.set_ylabel('Postion (m)',fontsize=14)
    ax.set_title('Uniform Circular Array Positions',fontsize=18)
    
    fig1,ax1 = plt.subplots()
    for ph in [0,15,30,45,60,75,90,105,120,135,150,165,180]:
        ph = int(ph)
        ax1.plot(np.rad2deg(theta),np.abs(lowTot[:,ph]),label=f'Phi={ph}')
    ax1.grid(True,which='both')
    ax1.legend(loc='upper right')
    ax1.set_title('Total')
    fig.set_size_inches(10,8)

    plt.show()