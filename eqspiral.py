import os, tempfile
import numpy as np

from CSXCAD  import ContinuousStructure
from openEMS import openEMS
from openEMS.physical_constants import *

def makeSpiral(FDTD,CSX,mesh,center,radii,alpha,h,hs,Nidx,excite):

    # Calculate top spiral
    r0 = radii[0]
    rmax = radii[1]
    N = np.log(rmax/r0)/(2*np.pi*alpha)
    a0 = np.pi/4
    a1 = np.linspace(a0,N*2*np.pi,np.round(100*N).astype(int))
    a2 = a1 + np.pi/2
    a3 = a2 + np.pi/2
    a4 = a3 + np.pi/2
    r = r0 * np.exp(alpha*a1)
    tp = np.zeros([2, 2*r.size+52])
    bp = np.zeros([2, 2*r.size+52])

    tp[0,:] = np.hstack((0,r,r[-1]*np.ones(50),np.flip(r),0))
    tp[1,:] = np.hstack((0,a1,np.linspace(a1[-1],a2[-1],50),np.flip(a2),0))
    bp[0,:] = np.hstack((0,r,r[-1]*np.ones(50),np.flip(r),0))
    bp[1,:] = np.hstack((0,a3,np.linspace(a3[-1],a4[-1],50),np.flip(a4),0))

    #create fixed lines for the simulation box and port
    mesh.AddLine('x', [center[0]-r0, center[0]-0.25, center[0]+0.25, center[0]+r0])
    mesh.AddLine('y', [center[1]-r0, center[1]-0.25, center[1]+0.25, center[1]+r0])
    bp_cart = np.zeros_like(tp)
    tp_cart = np.zeros_like(tp)
    tp_cart[0,:] = tp[0,:] * np.cos(tp[1,:])
    tp_cart[1,:] = tp[0,:] * np.sin(tp[1,:])
    bp_cart[0,:] = bp[0,:] * np.cos(bp[1,:])
    bp_cart[1,:] = bp[0,:] * np.sin(bp[1,:])
    # move coordinates to spiral center
    tp_cart[0,:] = center[0] + tp_cart[0,:]
    tp_cart[1,:] = center[1] + tp_cart[1,:]
    bp_cart[0,:] = center[0] + bp_cart[0,:]
    bp_cart[1,:] = center[1] + bp_cart[1,:]
    top = CSX.AddMetal(f'top_spiral{Nidx}')
    top.AddLinPoly(priority=10,norm_dir='z',elevation=h+hs,points=tp_cart,length=0)
    bot = CSX.AddMetal(f'bot_spiral{Nidx}')
    bot.AddLinPoly(priority=10,norm_dir='z',elevation=h-hs,points=bp_cart,length=0)

    substrate=CSX.AddMaterial(f'substrate{Nidx}')
    substrate.SetMaterialProperty(epsilon=4.2)
    start = [center[0]-rmax, center[1]-rmax, h-hs]
    stop =  [center[0]+rmax, center[1]+rmax, h+hs]
    substrate.AddBox(start=start, stop=stop, priority=0);


    # apply the excitation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    start=[center[0]-0.25, center[1]-0.25, h-hs]
    stop =[center[0]+0.25, center[1]+0.25, h+hs]
    port= FDTD.AddLumpedPort(priority=5,port_nr=Nidx,R=150,start=start,stop=stop,p_dir='z',excite=excite)
    return [CSX, FDTD, mesh, port]
    
def SimulateEmbeddedFarfield(freq,hs,h,centers,radii,theta,phi,eid=0):
    unit = 1e-3 # all length in mm
    f_start = 0.8 * freq
    f_stop = 1.2 *freq
    max_res = np.floor(C0 / (f_stop) / unit / 20) #cell size: lambda/20
    padding = max_res *20
    
    Sim_dir = os.path.join(tempfile.gettempdir(),'spiral_test')
    # size of the simulation box
    SimBox = np.array([padding+np.max(radii)*4+np.max(np.abs(centers[:,0])), padding+np.max(radii)*4+np.max(np.abs(centers[:,1])), padding+h+hs])


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
    mesh.AddLine('z', [0, h-hs, h+hs, SimBox[2]])

    # generate the spirals
    for idx in range(radii.shape[0]):
        this_r = radii[idx,:]
        this_center = centers[idx,:]
        if idx==0:
            excite=True
        else:
            excite=False
        [CSX,FDTD,mesh,port] = makeSpiral(FDTD,CSX,mesh,this_center,this_r,0.32,h,hs,idx,excite)
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
