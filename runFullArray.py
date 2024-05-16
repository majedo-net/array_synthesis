import os
#os.add_dll_directory('C:/openEMS/openEMS')
from antennaArray import AntennaArray
from simulation import Simulation
import numpy as np
from shutil import rmtree


if __name__ == '__main__':
    sims = []
    for idx in range(37):
        print('==================')
        print(f'Running Element {idx}')
        print('==================')
        ant_array = AntennaArray()
        freq = 6e9
        results_dir = '/results'
        cwd = os.getcwd()
        ant_array.generateRPSPositions(fmax=6e9,r=1,Nrps=3)
        ant_array.initDipoleElements(freq,orientation='z',excite_idx = [idx])
        theta = np.linspace(0, np.pi, 181)
        phi = np.linspace(0, 2*np.pi, 361)
        sim = Simulation(freq,theta,phi,ant_array,results_dir,id_=idx)
        sim.makeElementSims()
        sims.append(sim)
        sim.runSim()
        print(f'smn shape: {sim.smn[:,:,151].shape}\n')
        os.chdir(cwd)
        np.savetxt(f'{results_dir}/fullarray_sm{idx}.txt',sim.smn[:,:,151])
        np.savetxt(f'{results_dir}/full_s11_{idx}.txt',(sim.sfreqs,sim.s11[idx,:]))
        simdir = sim.simdir
        del sim
        del ant_array
        rmtree(simdir,ignore_errors=True)