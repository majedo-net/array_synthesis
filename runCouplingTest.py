import os
import numpy as np
import h5py
from multiprocessing import Pool
from shutil import rmtree
from antennaArray import AntennaArray
from simulation import Simulation

def run_simulation(Nr_idx):
    idx, fc = Nr_idx
    print('==================')
    print(f'Running Element fc = {fc}GHz')
    print('==================')
    
    ant_array = AntennaArray()
    freq = 6e9
    results_dir = f'/results/'
    cwd = os.getcwd()
    
    ant_array.generateRPSPositions(fmax=fc, r=1, Nrps=1)
    ant_array.excite_idx = idx
    ant_array.initDipoleElements(freq, orientation='y')
    
    theta = np.linspace(0, np.pi, 181)
    phi = np.linspace(0, 2*np.pi, 361)
    
    sim = Simulation(freq, theta, phi, ant_array, results_dir, id_=idx)
    sim.makeElementSims()
    sim.runSim()
    
    os.chdir(cwd)
    
    # Save HDF5 results
    with h5py.File(f'{results_dir}/fc{idx}.hdf5','w') as f:
        f.attrs['x'] = sim.array.elements[0].x
        f.attrs['y'] = sim.array.elements[0].y
        f.attrs['ant_type'] = sim.array.elements[0].ant_type
        f.attrs['fc'] = fc
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

# This is now the half wavelength freq in ghz
Nrs = [3,4,5,5.5,6,6.5]
Nr_idx_list = [(idx,Nr) for idx,Nr in enumerate(Nrs)]

# Create a Pool of workers to run simulations in parallel
with Pool(int(os.cpu_count()/4)) as pool:
    pool.map(run_simulation, Nr_idx_list)
