import pyswarms as ps
import random
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
    mins = np.ones(6)
    maxes = np.ones(6)
    bounds = (mins, maxes)
    mins[0],maxes[0] = 0.1,2 # r bounds
    mins[1],maxes[1] = -20,20 # rotation angle bounds
    mins[2],maxes[2] = -20,20 # rotation angle bounds
    mins[3],maxes[3] = -20,20 # rotation angle bounds
    mins[4],maxes[4] = -20,20 # rotation angle bounds
    mins[5],maxes[5] = 10,50 # Ground plane height bounds
    optimizer = ps.single.GlobalBestPSO(n_particles=n,dimensions=6, options=options, bounds=bounds)
    return optimizer

def compute_costs(pos):
    rs = pos[:,0]
    thetas = pos[:,1:4]
    hs = pos[:,5]
    costs = np.ones(len(rs))
    for n in range(len(rs)):
        hs[n] = int(hs[n])
        costs[n] = cost_function(rs[n],thetas[n],hs[n])
    return costs 

def check_completion(paths):
    not_complete = True
    while(not_complete):
        not_complete = False
        for p in paths:
            try:
                exists = os.path.exists(p)
                if exists==False:
                    not_complete=True
            except:
                pass

def check_element_patterns(spirads,h,freq,Np=0):
    if not os.path.exists('/results/ff'):
        os.mkdir('/results/ff')
    new_patterns = []

    # some spiral parameters. Load these from a file eventually...
    r0 = 5 #mm
    alpha = 0.32
    spirads = np.around(spirads*1000,0)
    freq = np.around(freq/1e6,-1)
    print('=============================================\n')
    print(f'New Simulations required for this iteration:\n')
    for rad in np.unique(spirads):
        if rad < 10:
            raise ArrayElementException(message='Element radius <10cm, 1000 cost assigned')
        patt_path = f'/results/ff/farfieldspiral_rad_{int(np.around(rad,0))}_freq_{int(freq)}_height_{h}.csv'
        # if the pattern file already exists, just copy it into the element pattern array
        # if file size gets unwieldy we may need to only store uniques for this too
        if os.path.exists(patt_path):
            continue
        else:
            print(patt_path)
            new_patterns.append(patt_path)
            with open('/results/commands.txt',mode='a+') as f, open('/results/commands_ff.txt',mode='a+') as f2:
                f.write(f'octave --silent spiral.m {int(freq)}e6 {r0} {alpha} {h} {rad} {Np} {int(freq)}\n')
                f2.write(f'octave --silent spiral_ff.m {int(freq)}e6 {r0} {alpha} {h} {rad} {Np} {int(freq)}\n')
                Np =Np+ 1

    print('=============================================\n')
    return new_patterns

def fetch_element_patterns(spirads,freq,h):
    patterns = []
    for i,rad in enumerate(spirads):
        rad = rad*1000
        patt_path = f'/results/ff/farfieldspiral_rad_{int(np.around(rad,0))}_freq_{int(np.around(freq/1e6,-1))}_height_{h}.csv'
        patt=np.genfromtxt(patt_path,delimiter=',',max_rows=181)
        print(f'{patt_path}    ==== Found and parsed \n')
        patt = patt / np.amax(patt)
        patterns.append(patt)
    return patterns

def cost_function(r,thetas,h,plots=False):
    N = 3
    fmax = 2.3e9
    fmin = 0.5e9
    d = rps(fmax,r,N)
    h=int(h)
    xs, ys = circ_positions(d)
    xs, ys = rotate(xs,ys,thetas,d)
    circs = [None]*len(xs)
    points = list(zip(xs.ravel(),ys.ravel()))
    spirads = KDTree(points).query(points,k=[2])[0].flatten()
    spirads = np.around((spirads/2)*0.96,2)

    freq = fmax * (1+np.sin(60*np.pi/180))
    try:
        np1=check_element_patterns(spirads,h,freq)
        np2=check_element_patterns(spirads,h,fmin,Np=len(np1)+1)
        new_patts = np1 + np2
        invoke_openems()
    except ArrayElementException as e:
        print(e)
        cost = 100000
        return cost
    element_pattern_max = fetch_element_patterns(spirads,freq,h)
    element_pattern_min = fetch_element_patterns(spirads,fmin,h)
    lam_max = 3e8/freq
    lam_min = 3e8/fmin
    k_max = 2*np.pi/lam_max
    k_min = 2*np.pi/lam_min
    theta = np.linspace(-np.pi/2, np.pi/2, 181)
    phi = np.linspace(0, np.pi, 181)
    ArrF_max, Tot_max = array_factor(xs,ys,k_max,element_pattern_max,theta,phi)
    ArrF_min, Tot_min = array_factor(xs,ys,k_min,element_pattern_min,theta,phi)
    des_bw = 35
    sll = -5
    meas_bw = Beamwidth(theta,ArrF_max)
    costmax,peaksmax = BeamCost(des_bw,meas_bw,theta,phi,Tot_max)
    meas_bw = Beamwidth(theta,ArrF_min)
    costmin,peaksmin = BeamCost(des_bw,meas_bw,theta,phi,Tot_min)
    cost = costmin + costmax

    if not os.path.exists('/results/plots'):
        os.mkdir('/results/plots')
    idstring = '%08x' % random.randrange(16**8)
    makeArrayPlot(xs,ys,d,spirads,idstring)
    with open(f'/results/{idstring}_array.txt',mode='w') as f: 
        f.write(f'h: {h} \n\n xs: {xs} \n\n ys: {ys} \n\n spirads: {spirads}')

    makePatternPlots(theta,phi,ArrF_max,Tot_max,cost,freq,idstring,save=True)
    makePatternPlots(theta,phi,ArrF_min,Tot_min,cost,fmin,idstring,save=True)
    return cost

def invoke_openems():
    if os.path.exists('/results/commands.txt'):
        # Run time stepping
        proc=subprocess.Popen("parallel -j5 < /results/commands.txt",shell=True)
        while proc.poll() is None:
            pass
        # Run far field calculation
        print('====================================\n')
        print('Starting Far Fields\n')
        print('====================================\n')
        proc=subprocess.Popen("parallel -j15 < /results/commands_ff.txt",shell=True)
        while proc.poll() is None:
            pass
        os.remove('/results/commands.txt')
        os.remove('/results/commands_ff.txt')
        return
    else:
        print('No new sims\n')
        return

if __name__ == '__main__':
    n_particles =3 
    n_iterations =3 
    swarm=init_swarm(n_particles)
    cost,pos = swarm.optimize(compute_costs,n_iterations)
    cost_function(pos[0],pos[1:-1],n_particles,plots=True)
    with open('results/optimized',mode='w') as f: f.write(f'cost: {cost} \n pos: {pos}')
