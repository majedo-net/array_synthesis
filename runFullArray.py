from antennaArray import AntennaArray
from simulation import Simulation
import numpy as np


if __name__ == '__main__':
    ant_array = AntennaArray()
    freq = 6e9
    ant_array.generateRPSPositions(fmax=6e9,r=1,Nrps=3)
    ant_array.initDipoleElements(freq,orientation='z')
    theta = np.linspace(0, np.pi, 181)
    phi = np.linspace(0, 2*np.pi, 361)
    sim = Simulation(freq,theta,phi,ant_array,'test')
    sim.makeElementSims()
    sim.runSim()
    print(f'smn shape: {sim.smn[:,:,151].shape}\n')
    np.savetxt('test/fullarray_smn.txt',sim.smn[:,:,151])
    print('break')