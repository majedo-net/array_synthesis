import os
import numpy as np
from multiprocessing import Pool
from shutil import rmtree
from antennaArray import AntennaArray
from simulation import Simulation

def run_simulation(Nr_idx):
    Nr, idx = Nr_idx
    print('==================')
    print(f'Running Element {idx} Tess Array')
    print('==================')
    
    ant_array = AntennaArray()
    freq = 6e9
    results_dir = f'/results/nr{Nr}'
    cwd = os.getcwd()
    
    ant_array.generateRPSPositions(fmax=6e9, r=1, Nrps=3)
    ant_array.excite_idx = idx
    ant_array.initDipoleElements(freq, orientation='z')
    tess_array = ant_array.getNearestNeighborSA(idx, Nr)
    
    theta = np.linspace(0, np.pi, 181)
    phi = np.linspace(0, 2*np.pi, 361)
    
    sim = Simulation(freq, theta, phi, tess_array, results_dir, id_=(Nr*100 +idx))
    sim.makeElementSims()
    sim.runSim()
    
    os.chdir(cwd)
    
    nids = [el.id for el in sim.array.elements]
    print(nids)
    
    np.savetxt(f'{results_dir}/sm{idx}.txt', (nids, sim.smn[:, 151]))
    np.savetxt(f'{results_dir}/s11_{idx}.txt', (sim.sfreqs, sim.s11))
    
    simdir = sim.simdir
    del sim
    del ant_array, tess_array
    rmtree(simdir, ignore_errors=True)

Nrs = [4,6,8,10,12,14,16,18,20,24,28,32,36]
Nr_idx_list = [(Nr, idx) for Nr in Nrs for idx in range(37)]

# Create a Pool of workers to run simulations in parallel
with Pool() as pool:
    pool.map(run_simulation, Nr_idx_list)
