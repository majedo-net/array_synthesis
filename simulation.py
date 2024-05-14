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
        self.rdir=results_dir_
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


        
