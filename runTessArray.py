import os
#os.add_dll_directory('C:/openEMS/openEMS')
from antennaArray import AntennaArray
from simulation import Simulation
import numpy as np
from shutil import rmtree


if __name__ == '__main__':
    Nrs = [4,6,8,10,12,14,16]
    for Nr in Nrs:
        for idx in range(37):
            print('==================')
            print(f'Running Element {idx} Tess Array')
            print('==================')
            ant_array = AntennaArray()
            freq = 6e9
            results_dir = f'/results/nr{Nr}'
            cwd = os.getcwd()
            ant_array.generateRPSPositions(fmax=6e9,r=1,Nrps=3)
            ant_array.excite_idx = idx
            ant_array.initDipoleElements(freq,orientation='y')
            tess_array = ant_array.getNearestNeighborSA(idx,Nr)
            theta = np.linspace(0, np.pi, 181)
            phi = np.linspace(0, 2*np.pi, 361)
            sim = Simulation(freq,theta,phi,tess_array,results_dir,id_=idx)
            sim.makeElementSims()
            sim.runSim()
            os.chdir(cwd)
            nids = []
            for el in sim.array.elements: nids.append(el.id)
            print(nids)
            np.savetxt(f'{results_dir}/sm{idx}.txt',(nids,sim.smn[:,151]))
            np.savetxt(f'{results_dir}/s11_{idx}.txt',(sim.sfreqs,sim.s11))
            simdir = sim.simdir
            del sim
            del ant_array, tess_array
            rmtree(simdir,ignore_errors=True)