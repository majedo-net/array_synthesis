from antenna import Antenna
import os 
os.add_dll_directory('C:/openEMS/openEMS')
from CSXCAD import ContinuousStructure
from openEMS import openEMS
from openEMS.physical_constants import *

class PatchAntenna(Antenna):
    def __init__(self,id_,L_,W_,x_,y_,hs_,epsr_):
        super().__init__(id_,'Patch')
        self.L = L_
        self.W = W_
        self.x = x_
        self.y = y_
        self.hs = hs_
        self.epsr=epsr_

    def makeSim(self,FDTD,CSX,mesh,excite,max_res):
        patch = CSX.AddMetal(f'patch{self.id}')
        start = [self.x-self.W/2, self.y-self.L/2, self.hs]
        stop = [self.x+self.W/2, self.y+self.L/2, self.hs]
        patch.AddBox(priority=10, start=start, stop=stop)
        mesh.AddLine('x',[start[0],stop[0]])
        mesh.AddLine('y',[start[1],stop[1]])
        FDTD.AddEdges2Grid(dirs='xy',properties=patch,metal_edge_res=max_res)

        substrate=CSX.AddMaterial(f'substrate{self.id}')
        substrate.SetMaterialProperty(epsilon=self.epsr)
        start = [self.x-self.W, self.y-self.L, 0]
        stop = [self.x+self.W, self.y+self.L, self.hs]
        substrate.AddBox(start=start, stop=stop, priority=0)

        ground =CSX.AddMetal(f'ground{self.id}')
        start = [self.x-self.W, self.y-self.L, 0]
        stop = [self.x+self.W, self.y+self.L, 0]
        ground.AddBox(start=start,stop=stop,priority=10)
        FDTD.AddEdges2Grid(dirs='xy',properties=ground)

        mesh.AddLine('z', np.linspace(0,self.hs,5))

        # apply the excitation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        start=[self.x-0.2,self.y-0.22*self.L,0]
        stop =[self.x+0.2,self.y-0.18*self.L,self.hs]
        port= FDTD.AddLumpedPort(priority=5,port_nr=self.id,R=50,start=start,stop=stop,p_dir='z',excite=self.excite)

        mesh.AddLine('x',[start[0],stop[0]])
        mesh.AddLine('y',[start[1],stop[1]])
        return [CSX, FDTD, mesh, port]

