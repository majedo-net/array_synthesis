import os, tempfile
import numpy as np

from CSXCAD  import ContinuousStructure
from openEMS import openEMS
from openEMS.physical_constants import *

def makePatch(FDTD,CSX,mesh,center,L,W,hs,max_res,Nidx,excite):

    patch = CSX.AddMetal(f'patch{Nidx}')
    start = [center[0]-W/2, center[1]-L/2, hs]
    stop = [center[0]+W/2, center[1]+L/2, hs]
    patch.AddBox(priority=10, start=start,stop=stop)
    FDTD.AddEdges2Grid(dirs='xy',properties=patch,metal_edge_res=max_res)

    substrate=CSX.AddMaterial(f'substrate{Nidx}')
    substrate.SetMaterialProperty(epsilon=4.2)
    start = [center[0]-W, center[1]-L, 0]
    stop = [center[0]+W, center[1]+L, hs]
    substrate.AddBox(start=start, stop=stop, priority=0)

    ground =CSX.AddMetal(f'ground{Nidx}')
    start = [center[0]-W, center[1]-L, 0]
    stop = [center[0]+W, center[1]+L, 0]
    ground.AddBox(start=start,stop=stop,priority=10)
    FDTD.AddEdges2Grid(dirs='xy',properties=ground)

    mesh.AddLine('z', np.linspace(0,hs,5))

    # apply the excitation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    start=[center[0]-1,center[1]-0.25*L,0]
    stop =[center[0]+1,center[1]-0.25*L,hs]
    port= FDTD.AddLumpedPort(priority=5,port_nr=Nidx,R=50,start=start,stop=stop,p_dir='z',excite=excite)
    return [CSX, FDTD, mesh, port]
    
def SimulateEmbeddedFarfield(freq,hs,centers,L,W,theta,phi,eid=0):
    unit = 1e-3 # all length in mm
    f_start = 0.8 * freq
    f_stop = 1.2 *freq
    max_res = np.floor(C0 / (f_stop) / unit / 20) #cell size: lambda/20
    padding = max_res *20
    
    Sim_dir = os.path.join(tempfile.gettempdir(),'spiral_test')
    # size of the simulation box
    SimBox = np.array([padding+L*2+np.max(np.abs(centers[:,0])), padding+L+np.max(np.abs(centers[:,1])), padding+hs])


    ## setup FDTD parameter & excitation function
    FDTD = openEMS(EndCriteria=1e-4,NrTS=50000)
    FDTD.SetGaussExcite(0.5*(f_start+f_stop),0.5*(f_stop-f_start))
    FDTD.SetBoundaryCond(['PML_8', 'PML_8', 'PML_8', 'PML_8', 'PEC', 'PML_8']) # boundary conditions

    CSX = ContinuousStructure()
    FDTD.SetCSX(CSX)
    mesh = CSX.GetGrid()
    mesh.SetDeltaUnit(unit)
    ports = []
    
    # create fixed lines for the simulation box and port
    mesh.AddLine('x', [-SimBox[0]/2, SimBox[0]/2])
    mesh.AddLine('y', [-SimBox[1]/2, SimBox[1]/2])
    mesh.AddLine('z', [0, hs, SimBox[2]])

    # generate the spirals
    for idx in range(centers.shape[0]):
        this_center = centers[idx,:]
        if idx==0:
            excite=True
        else:
            excite=False
        [CSX,FDTD,mesh,port] = makePatch(FDTD,CSX,mesh,this_center,L,W,hs,max_res,idx,excite)
        ports.append(port)

    mesh.SmoothMeshLines('all',max_res,1.4)
    nf2ff = FDTD.CreateNF2FFBox()
    CSX.Write2XML(f'/results/csx{eid}.xml')
    FDTD.Run(Sim_dir, cleanup=True)
    ffres = nf2ff.CalcNF2FF(Sim_dir,freq,theta,phi,center=centers[0,:])
    sfreqs = np.linspace(f_start,f_stop,101)
    ports[0].CalcPort(Sim_dir,sfreqs)
    s11 = ports[0].uf_ref / ports[0].uf_inc
    s11 = np.abs(s11)
    s11_db = 20.0*np.log10(s11)
    s11f = s11[50]
    E_norm = ffres.E_norm[0]/np.max(ffres.E_norm[0]) 

    return E_norm,s11f,s11_db,sfreqs


if __name__ == '__main__':
    centers = np.array([[0, 0], [200, 100],[100, -100],[-100, -100],[-100, 100]])
    radii = np.array([[5,50],[5,40],[5,30],[5,20],[5,35]])
    freq = 4e9
    hs = 0.2 # substrate thickness
    h = 10 # cavity height
    theta = np.linspace(0, np.pi, 181)
    phi = np.linspace(0, 2*np.pi, 361)

    En,s11db,sfreqs = SimulateEmbeddedFarfield(freq,hs,h,centers,radii,theta,phi)

    print(En.shape)
