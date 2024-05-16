from antenna import Antenna
from CSXCAD import ContinuousStructure
from openEMS import openEMS
from openEMS.physical_constants import *

class DipoleAntenna(Antenna):
    def __init__(self,id_,f0,x_,y_,orientation_='x',excite_=True):
        super().__init__(id_,'Dipole',excite_=excite_)
        lamb = 3e8/f0
        self.L = 0.445*lamb*1000
        self.x= x_
        self.y= y_
        self.orientation = orientation_

    def makeSim(self,FDTD,CSX,mesh,excite,max_res):
        dipole = CSX.AddMetal(f'dipole{self.id}')
        FDTD.SetBoundaryCond(['PML_8', 'PML_8', 'PML_8', 'PML_8', 'PML_8', 'PML_8'])
        if self.orientation == 'x':
            start = [self.x-self.L/2, self.y-0.5, 3-0.5]
            stop = [self.x+self.L/2, self.y+0.5, 3+0.5]
            pstart = [self.x-max_res/2, self.y-0.5, 3-0.5]
            pstop = [self.x+max_res/2, self.y+0.5, 3+0.5]
            mesh.AddLine('x',[start[0],start[0]+(stop[0]-start[0])/3,start[0]+2*(stop[0]-start[0])/3,stop[0],pstart[0],pstop[0]])
            mesh.AddLine('y',[start[1]])
            mesh.AddLine('z',[start[2]])
        if self.orientation == 'y':
            start = [self.x-0.5, self.y-self.L/2, 3-0.5]
            stop = [self.x+0.5, self.y+self.L/2, 3+0.5]
            pstart= [self.x-0.5, self.y-max_res/2, 3-0.5]
            pstop = [self.x+0.5, self.y+max_res/2, 3+0.5]
            mesh.AddLine('y',[start[1],stop[1],start[1]+(stop[1]-start[1])/3,start[1]+2*(stop[1]-start[1])/3,pstart[1],pstop[1]])
            mesh.AddLine('x',[start[0]])
            mesh.AddLine('z',[start[2]])
        if self.orientation == 'z':
            start = [self.x-0.5, self.y-0.5, -self.L/2]
            stop = [self.x+0.5,self.y+0.5, self.L/2]
            pstart = [self.x-0.5, self.y-0.5, -max_res/2]
            pstop = [self.x+0.5, self.y+0.5, max_res/2]
            mesh.AddLine('z',[start[2],start[2]+(stop[2]-start[2])/3,start[2]+2*(stop[2]-start[2])/3,stop[2],pstart[2],pstart[2]])
            mesh.AddLine('x',[start[0]])
            mesh.AddLine('y',[start[1]])
        dipole.AddBox(priority=5,start=start,stop=stop)
        dirs = 'xyz'.strip(self.orientation)
        FDTD.AddEdges2Grid(dirs='all',properties=dipole)
        port = FDTD.AddLumpedPort(priority=1,port_nr=self.id,R=70,start=pstart,stop=pstop,p_dir=self.orientation,excite=self.excite)
        return [CSX, FDTD, mesh, port]
        
        



