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
    mesh.AddLine('x',[start[0],stop[0]])
    mesh.AddLine('y',[start[1],stop[1]])
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

    start=[center[0]-0.2,center[1]-0.22*L,0]
    stop =[center[0]+0.2,center[1]-0.18*L,hs]
    port= FDTD.AddLumpedPort(priority=5,port_nr=Nidx,R=50,start=start,stop=stop,p_dir='z',excite=excite)

    mesh.AddLine('x',[start[0],stop[0]])
    mesh.AddLine('y',[start[1],stop[1]])
    return [CSX, FDTD, mesh, port]
    
def SimulateEmbeddedFarfield(freq,hs,centers,L,W,theta,phi,eid=0):
    unit = 1e-3 # all length in mm
    f_start = 0.8 * freq
    f_stop = 1.2 *freq
    max_res = np.floor(C0 / (f_stop) / unit / 20) #cell size: lambda/20
    padding = max_res *20
    
    Sim_dir = os.path.join(tempfile.gettempdir(),f'patch_sim{eid}')
    # size of the simulation box
    SimBox = np.array([padding+W*2+np.max(np.abs(centers[:,0])), padding+L*2+np.max(np.abs(centers[:,1])), padding+hs])


    ## setup FDTD parameter & excitation function
    FDTD = openEMS(EndCriteria=1e-4,NrTS=500000)
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

    # generate the patches
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
    ffres = nf2ff.CalcNF2FF(Sim_dir,freq,theta,phi)
    sfreqs = np.linspace(f_start,f_stop,301)
    ports[0].CalcPort(Sim_dir,sfreqs)
    s11 = ports[0].uf_ref / ports[0].uf_inc
    zin = ports[0].uf_tot / ports[0].if_tot
    s11 = np.abs(s11)
    s11_db = 20.0*np.log10(s11)
    s11f = s11[50]
    E_norm = ffres.E_norm[0]/np.max(ffres.E_norm[0]) 

    return E_norm,s11f,s11_db,sfreqs, zin


if __name__ == '__main__':
    centers = np.array([[0, 0]])
    L = 10.2
    W = 15.5
    hs = 1.5
    freqs = [6e9]
    hs = 1.5 # substrate thickness
    theta = np.linspace(0, np.pi, 181)
    phi = np.linspace(0, 2*np.pi, 361)

    En,s11f,s11db,sfreqs,zin = SimulateEmbeddedFarfield(freqs[0],hs,centers,L,W,theta,phi)

    print(zin.shape)

    np.savetxt('/results/s11_single_patch.txt',(sfreqs,s11db,np.real(zin),np.imag(zin)))
    np.savetxt('/results/ff_single_patch.txt',En)