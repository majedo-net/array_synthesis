from antennaArray import AntennaArray
from simulation import Simulation
import numpy as np


if __name__ == '__main__':
    for idx in range(37):
        print('==================')
        print(f'Running Element {idx}')
        print('==================')
        ant_array = AntennaArray()
        freq = 6e9
        ant_array.generateRPSPositions(fmax=6e9,r=1,Nrps=3)
        ant_array.initDipoleElements(freq,orientation='z',excite_idx = [idx])
        theta = np.linspace(0, np.pi, 181)
        phi = np.linspace(0, 2*np.pi, 361)
        sim = Simulation(freq,theta,phi,ant_array,'/results',id_=idx)
        sim.makeElementSims()
        sim.runSim()
        print(f'smn shape: {sim.smn[:,:,151].shape}\n')
        np.savetxt(f'/results/fullarray_sm{idx}.txt',sim.smn[:,:,151])
        np.savetxt(f'/results/full_s11_{idx}.txt',(sim.sfreqs,sim.s11[idx,:]))
        del sim
        del ant_array