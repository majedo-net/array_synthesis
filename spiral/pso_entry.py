import pyswarms as ps
import numpy as np
import subprocess,os,fileinput,re,shutil
regex = re.compile(r'\d+')

def init_swarm(n):
    options = {'c1': 0.5, 'c2': 0.3, 'w':0.3}
    mins = np.ones(3)
    maxes = np.ones(3)
    bounds = (mins, maxes)
    mins[0],maxes[0] = 2,10 # r0 bounds
    mins[1],maxes[1] = 0.1,0.8 # alpha bounds
    mins[2],maxes[2] = 31,40 # h bounds
    optimizer = ps.single.GlobalBestPSO(n_particles=n,dimensions=3, options=options, bounds=bounds)
    return optimizer

def cost_function(pos,rmax):
    r0 = pos[:,0]
    alpha = pos[:,1]
    h = pos[:,2]
    freqs = [
            [0.3e9,0.6e9,'a'],
            [0.6e9,1.5e9,'b'],
            [1.5e9,2.5e9,'c']]
    if os.path.isfile('commands.txt'): os.remove('commands.txt')
    with open('commands.txt',mode='a') as commands:    
        for p in range(len(r0)):
            for f in freqs:
                   commands.write(f'octave --silent spiral.m {f[0]} {f[1]} {r0[p]} {alpha[p]} {h[p]} {rmax} {p} {f[2]}\n') 
    shutil.copyfile('commands.txt','results/commands.txt')
    invoke_openems()
    # Concatenate record keeping file
    cost = np.ones(len(r0))
    for n in range(len(r0)):
        inputfiles=[f"results/{x}" for x in os.listdir('results') if re.match(f"GZresults{n}[a-c]",x)]
        with open(f"results/GZBresults{n}.csv",'w') as fout, fileinput.input(inputfiles) as fin:
            for line in fin: fout.write(line)
        for f in inputfiles:os.remove(f)
        d = np.loadtxt(f'results/GZBresults{n}.csv',delimiter=',')
        zreal_std = np.std(d[:,1])
        zimag_std = np.std(d[:,2])
        g_std = np.std(d[:,3])
        s11_max = np.max(d[:,4])
        s11_ave = np.mean(d[:,4])
        cost[n] = s11_max + s11_ave
    return cost

def invoke_openems():
   proc=subprocess.Popen("parallel -j15 < commands.txt",shell=True)
   while proc.poll() is None:
       pass

   return

if __name__ == '__main__':
    n_particles = 7
    n_iterations = 10
    swarm=init_swarm(n_particles)
    cost,pos = swarm.optimize(cost_function,n_iterations,rmax=250)
    with open('results/optimized') as f: f.write(f'cost: {cost} \n pos: {pos}')
