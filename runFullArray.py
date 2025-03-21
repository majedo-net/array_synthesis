import os
#os.add_dll_directory('C:/openEMS/openEMS')
from antennaArray import AntennaArray
from simulation import Simulation
import h5py
import numpy as np
from shutil import rmtree


if __name__ == '__main__':
    for idx in range(37):
        print('==================')
        print(f'Running Element {idx}')
        print('==================')
        ant_array = AntennaArray()
        freq = 6e9
        results_dir = '/results'
        cwd = os.getcwd()
        ant_array.generateRPSPositions(fmax=6e9,r=1,Nrps=3)
        ant_array.excite_idx = idx
        #ant_array.initPatchElements()
        #ant_array.initDipoleElements(freq,orientation='y')
        ant_array.initSpiralElements(r0=1,rmax=10,h=10)
        theta = np.linspace(0, np.pi, 181)
        phi = np.linspace(0, 2*np.pi, 361)
        sim = Simulation(freq,theta,phi,ant_array,results_dir,id_=idx)
        sim.makeElementSims()
        sim.runSim()
        os.chdir(cwd)
        with h5py.File(f'{results_dir}/element{idx}.hdf5','w') as f:
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
        del ant_array
        rmtree(simdir,ignore_errors=True)
