import os
import numpy as np
import h5py
from multiprocessing import Pool
from shutil import rmtree
from antennaArray import AntennaArray
from simulation import Simulation

def run_simulation(Nr_idx):
    Nr, idx = Nr_idx
    print('==================')
    print(f'Running Element {idx} Multithread Tess Array')
    print('==================')
    
    ant_array = AntennaArray()
    freq = 6e9
    results_dir = f'/results/nr{Nr}'
    cwd = os.getcwd()
    
    ant_array.generateRPSPositions(fmax=6e9, r=1, Nrps=3)
    ant_array.excite_idx = idx
    ant_array.initDipoleElements(freq, orientation='z')
    #ant_array.initSpiralElements(r0=1,rmax=10,h=10)
    tess_array = ant_array.getNearestNeighborSA(idx, Nr)
    
    theta = np.linspace(0, np.pi, 181)
    phi = np.linspace(0, 2*np.pi, 361)
    
    sim = Simulation(freq, theta, phi, tess_array, results_dir, id_=(Nr*100 +idx))
    sim.makeElementSims()
    sim.runSim()
    
    os.chdir(cwd)
    
    nids = [el.id for el in sim.array.elements]
    print(f'{idx}: {nids}')
    
    # Save HDF5 results
    with h5py.File(f'{results_dir}/element{idx}.hdf5','w') as f:
        nids = f.create_dataset('nids',data=nids)
        f.attrs['x'] = sim.array.elements[0].x
        f.attrs['y'] = sim.array.elements[0].y
        f.attrs['ant_type'] = sim.array.elements[0].ant_type
        smn = f.create_dataset('smn',data=sim.smn,compression='gzip')
        sfreqs = f.create_dataset('sfreqs',data=sim.sfreqs,compression='gzip')
        eth = f.create_dataset('eth',data=sim.eth,compression='gzip')
        eth.attrs['freq'] = sim.freq
        eph = f.create_dataset('eph',data=sim.eph,compression='gzip')
        eph.attrs['freq'] = sim.freq
        dmax = f.create_dataset('dmax',data=sim.dmax,compression='gzip')
        dmax.attrs['freq'] = sim.freq
    
    simdir = sim.simdir
    del sim
    del ant_array, tess_array
    rmtree(simdir, ignore_errors=True)
    print('returned')

Nrs1 = [[2],[4],[6],[8],[10],[12],[14]]
#Nrs1 = [[12],[14],[16],[18],[20],[24],[28],[32],[36]]
for Nrs in Nrs1:
    print(f'Starting Nrs = {Nrs}')
    Nr_idx_list = [(Nr, idx) for Nr in Nrs for idx in range(37)]
    print(Nr_idx_list)
    # Create a Pool of workers to run simulations in parallel
    with Pool(int(os.cpu_count()/4)) as pool:
        pool.map_async(run_simulation, Nr_idx_list)
        pool.close()
        pool.join()
