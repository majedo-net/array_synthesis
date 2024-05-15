import numpy as np
from CSXCAD import ContinuousStructure
from openEMS import openEMS
from openEMS.physical_constants import *
import os,tempfile

class Simulation():
    '''
    class to hold data and results defining a simulation
    '''

    def __init__(self,freq_,thetas_,phis_,array_,results_dir_,id_=0):
        self.unit=1e-3
        self.id = id_
        self.freq=freq_
        self.thetas=thetas_
        self.phis=phis_
        self.array=array_
        self.simdir=os.path.join(tempfile.gettempdir(),f'patch_sim{id_}')
        self.results_dir=results_dir_
        self.f_start=0.8*self.freq
        self.f_stop=1.2*self.freq
        self.max_res=np.floor(C0 / self.f_stop / self.unit / 20)
        self.padding = self.max_res * 20

        self.FDTD = openEMS(EndCriteria=1e-4,NrTS=5e6)
        self.FDTD.SetGaussExcite(0.5*(self.f_start+self.f_stop),0.5*(self.f_stop-self.f_start))
        self.FDTD.SetBoundaryCond(['PML_8', 'PML_8', 'PML_8', 'PML_8', 'PEC', 'PML_8']) # boundary conditions
        self.SimBox = np.array([self.padding+self.array.xmax, self.padding+self.array.ymax, self.padding+self.array.zmax])
        self.CSX = ContinuousStructure()
        self.FDTD.SetCSX(self.CSX)
        self.mesh = self.CSX.GetGrid()
        self.mesh.SetDeltaUnit(self.unit)
        self.ports = []

        self.mesh.AddLine('x', [-self.SimBox[0]/2, self.SimBox[0]/2])
        self.mesh.AddLine('y', [-self.SimBox[1]/2, self.SimBox[1]/2])
        self.mesh.AddLine('z', [0,self.SimBox[2]])

    def makeElementSims(self):
        for el in self.array.elements:
            [CSX, FDTD, mesh, port] = el.makeSim(self.FDTD,self.CSX,self.mesh,el.excite,self.max_res)
            self.CSX = CSX
            self.FDTD = FDTD
            self.mesh = mesh
            self.ports.append(port)

        if self.array.elements[0].ant_type=='Dipole': 
            self.mesh.AddLine('z', [-self.SimBox[2]/2,self.SimBox[2]/2])

        self.mesh.SmoothMeshLines('all',self.max_res,1.4)

    def runSim(self):
        self.nf2ff = self.FDTD.CreateNF2FFBox()
        self.CSX.Write2XML(f'{self.results_dir}/csx{self.id}.xml')
        self.FDTD.Run(self.simdir, cleanup=True)
        self.ffres = self.nf2ff.CalcNF2FF(self.simdir,self.freq,self.thetas,self.phis)
        sfreqs = np.linspace(self.f_start,self.f_stop,301)
        self.smn = np.zeros([len(self.ports),len(self.ports),sfreqs.shape[0]],dtype=np.complex128)
        for idx in range(len(self.ports)):
            print(f'calculating port {idx}')
            self.ports[idx].CalcPort(self.simdir,sfreqs)
    
        for midx in range(len(self.ports)):
            for nidx in range(len(self.ports)):
                pm = self.ports[midx].uf_ref
                pn = self.ports[nidx].uf_inc
                self.smn[midx,nidx,:] = pm/pn





        
