import pyswarms as ps
import numpy as np
import subprocess,os,fileinput,re,shutil
from hex_rps import *
from array_funcs import *
regex = re.compile(r'\d+')

def init_swarm(n):
    options = {'c1': 0.5, 'c2': 0.3, 'w':0.3}
    mins = np.ones(5)
    maxes = np.ones(5)
    bounds = (mins, maxes)
    mins[0],maxes[0] = 0.1,1.3 # r bounds
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
        costs[n] = cost_function(rs[n],thetas[n])
    return costs 

def cost_function(r,thetas):
    N = 3
    fmax = 3e9
    fmin = 0.35e9
    d = rps(fmax,r,N)
    xs, ys = hex_positions(d)
    xs, ys = rotate(xs,ys,thetas,d)
    element_pattern_max = np.ones(len(xs))
    element_pattern_min = np.ones(len(xs))
    freq = fmax * (1+np.sin(60*np.pi/180))
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

   return

if __name__ == '__main__':
    n_particles = 20
    n_iterations = 200
    swarm=init_swarm(n_particles)
    cost,pos = swarm.optimize(compute_costs,n_iterations)
    with open('results/optimized',mode='w') as f: f.write(f'cost: {cost} \n pos: {pos}')
