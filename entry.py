import pyswarms as ps
import numpy as np
from matplotlib.pyplot import *
from matplotlib.patches import Circle
from scipy.spatial import KDTree
import scipy
import subprocess,os,fileinput,re,shutil
from hex_rps import *
from circ_rps import *
from array_funcs import *
from ArrayElementException import ArrayElementException
regex = re.compile(r'\d+')

def init_swarm(n):
    options = {'c1': 0.5, 'c2': 0.3, 'w':0.3}
    mins = np.ones(5)
    maxes = np.ones(5)
    bounds = (mins, maxes)
    mins[0],maxes[0] = 0.1,2 # r bounds
    mins[1],maxes[1] = -20,20 # rotation angle bounds
    mins[2],maxes[2] = -20,20 # rotation angle bounds
    mins[3],maxes[3] = -20,20 # rotation angle bounds
    mins[4],maxes[4] = -20,20 # rotation angle bounds
    optimizer = ps.single.GlobalBestPSO(n_particles=n,dimensions=5, options=options, bounds=bounds)
    return optimizer

def compute_costs(pos):
    rs = pos[:,0]
    thetas = pos[:,1:4]
    costs = np.ones(len(rs))
    for n in range(len(rs)):
        costs[n] = cost_function(rs[n],thetas[n],n)
    return costs 

def check_element_patterns(spirads,freq):
    if not os.path.exists('/results/ff'):
        os.mkdir('/results/ff')
    patterns = []
    Np = 0 # an index for keeping track of any new simulations we have to run

    # some spiral parameters. Load these from a file eventually...
    r0 = 5 #mm
    alpha = 0.32
    h = 30 #mm
    spirads = np.around(spirads*1000,0)
    freq = np.around(freq/1e6,-1)
    print('=============================================\n')
    print(f'New Simulations required for this iteration:\n')
    for rad in np.unique(spirads):
        if rad < 10:
            raise ArrayElementException(message='Element radius <10cm, 1000 cost assigned')
        patt_path = f'/results/ff/farfieldspiral_rad_{int(np.around(rad,0))}_freq_{int(freq)}.csv'
        # if the pattern file already exists, just copy it into the element pattern array
        # if file size gets unwieldy we may need to only store uniques for this too
        if os.path.exists(patt_path):
            continue
        else:
            print(patt_path)
            with open('commands.txt',mode='a+') as f:
                f.write(f'octave --silent spiral.m {int(freq)}e6 {r0} {alpha} {h} {rad} {Np} \n')
                Np += 1

    print('=============================================\n')

def fetch_element_patterns(spirads,freq):
    patterns = []
    for i,rad in enumerate(spirads):
        rad = rad*1000
        patt_path = f'/results/ff/farfieldspiral_rad_{int(np.around(rad,0))}_freq_{int(np.around(freq/1e6,-1))}.csv'
        patterns.append(np.genfromtxt(patt_path,delimiter=','))
    return patterns

def cost_function(r,thetas,n_particle,plots=False):
    N = 3
    fmax = 1e9
    fmin = 0.5e9
    d = rps(fmax,r,N)
    xs, ys = circ_positions(d)
    xs, ys = rotate(xs,ys,thetas,d)
    circs = [None]*len(xs)
    points = list(zip(xs.ravel(),ys.ravel()))
    spirads = KDTree(points).query(points,k=[2])[0].flatten()
    spirads = np.around(spirads/2,2)
    if plots:
        fig,ax = subplots()
        ax.set_xlim(xmin=-np.max(d)*1.5,xmax=np.max(d)*1.5)
        xmax = 1.1*np.max(xs)+np.max(spirads)
        ymax = 1.1*np.max(ys)+np.max(spirads)
        ax.axis([-xmax, xmax, -ymax, ymax])
        dind = 1
        for i in range(len(xs)):
            circs[i] = Circle((xs[i],ys[i]),spirads[i],facecolor='blue',edgecolor='k')
            ax.add_patch(circs[i])
        ax.grid()
        ax.set_aspect('equal','box')
        fig.savefig('circles.png')
        print(np.unique(spirads))

    freq = fmax * (1+np.sin(60*np.pi/180))
    try:
        check_element_patterns(spirads,freq)
        invoke_openems()
        check_element_patterns(spirads,fmin)
        invoke_openems()
    except ArrayElementException as e:
        print(e)
        cost = 100000
        return cost
    element_pattern_max = fetch_element_patterns(spirads,freq)
    element_pattern_min = fetch_element_patterns(spirads,fmin)
    lam_max = 3e8/freq
    lam_min = 3e8/fmin
    k_max = 2*np.pi/lam_max
    k_min = 2*np.pi/lam_min
    theta = np.linspace(-np.pi/2, np.pi/2, 181)
    phi = np.linspace(0, np.pi, 181)
    ArrF_max = array_factor(xs,ys,k_max,element_pattern_max,theta,phi)
    ArrF_min = array_factor(xs,ys,k_min,element_pattern_min,theta,phi)
    des_bw = 35
    sll = -5
    meas_bw = Beamwidth(theta,ArrF_max)
    costmax = BeamCost(des_bw,meas_bw,sll,theta,phi,ArrF_max)
    meas_bw = Beamwidth(theta,ArrF_min)
    costmin = BeamCost(des_bw,meas_bw,sll,theta,phi,ArrF_min)
    cost = costmin + costmax

    return cost

def invoke_openems():
    proc=subprocess.Popen("parallel -j15 < commands.txt",shell=True)
    while proc.poll() is None:
        pass

    os.remove('commands.txt')
    return

if __name__ == '__main__':
    n_particles = 1
    n_iterations = 2
    swarm=init_swarm(n_particles)
    cost,pos = swarm.optimize(compute_costs,n_iterations)
    cost_function(pos[0],pos[1:-1],plots=True)
    #with open('results/optimized',mode='w') as f: f.write(f'cost: {cost} \n pos: {pos}')
