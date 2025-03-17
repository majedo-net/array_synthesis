from antenna import Antenna
import os 
import numpy as np
from CSXCAD import ContinuousStructure
from openEMS import openEMS
from openEMS.physical_constants import *

class SpiralAntenna(Antenna):
    def __init__(self,id_,x_,y_,r0_,rmax_,alpha_,h_,hs_,epsr_):
        super().__init__(id_,'Spiral')
        self.x = x_
        self.y = y_
        self.r0=r0_
        self.rmax=rmax_
        self.h=h_
        self.hs = hs_

        self.alpha=alpha_
        self.epsr=epsr_

    def makeSim(self,FDTD,CSX,mesh,excite,max_res):
        #calculate top spiral
        N = np.log(self.rmax/self.r0)/(2*np.pi*self.alpha)
        a0 = np.pi/4
        a1 = np.linspace(a0,N*2*np.pi,np.round(500*N).astype(int))
        a2 = a1 + np.pi/2
        a3 = a2 + np.pi/2
        a4 = a3 + np.pi/2
        r = self.r0 * np.exp(self.alpha*a1)
        tp = np.zeros([2, 2*r.size+52])
        bp = np.zeros([2, 2*r.size+52])

        
        tp[0,:] = np.hstack((0,r,r[-1]*np.ones(50),np.flip(r),0))
        tp[1,:] = np.hstack((0,a1,np.linspace(a1[-1],a2[-1],50),np.flip(a2),0))
        bp[0,:] = np.hstack((0,r,r[-1]*np.ones(50),np.flip(r),0))
        bp[1,:] = np.hstack((0,a3,np.linspace(a3[-1],a4[-1],50),np.flip(a4),0))

        h = self.h
        hs = self.hs
        mesh.AddLine('x', [self.x-self.r0, self.x-0.25, self.x+0.25, self.x+self.r0])
        mesh.AddLine('y', [self.y-self.r0, self.y-0.25, self.y+0.25, self.y+self.r0])
        mesh.AddLine('z', [h-hs, (h-hs)+2*hs/5, (h-hs)+3*hs/5, (h-hs)+4*hs/5, h+hs])
        
        bp_cart = np.zeros_like(tp)
        tp_cart = np.zeros_like(tp)
        tp_cart[0,:] = tp[0,:] * np.cos(tp[1,:])
        tp_cart[1,:] = tp[0,:] * np.sin(tp[1,:])
        bp_cart[0,:] = bp[0,:] * np.cos(bp[1,:])
        bp_cart[1,:] = bp[0,:] * np.sin(bp[1,:])
        # move coordinates to spiral center
        tp_cart[0,:] = self.x + tp_cart[0,:]
        tp_cart[1,:] = self.y + tp_cart[1,:]
        bp_cart[0,:] = self.x + bp_cart[0,:]
        bp_cart[1,:] = self.y + bp_cart[1,:]
        top = CSX.AddMetal(f'top_spiral{self.id}')
        top.AddLinPoly(priority=10,norm_dir='z',elevation=self.h+self.hs,points=tp_cart,length=0)
        bot = CSX.AddMetal(f'bot_spiral{self.id}')
        bot.AddLinPoly(priority=10,norm_dir='z',elevation=self.h-self.hs,points=bp_cart,length=0)

        substrate=CSX.AddMaterial(f'substrate{self.id}')
        substrate.SetMaterialProperty(epsilon=self.epsr)
        start = [self.x-self.rmax, self.y-self.rmax, self.h-self.hs]
        stop =  [self.x+self.rmax, self.y+self.rmax, self.h+self.hs]
        substrate.AddBox(start=start, stop=stop, priority=0);
        mesh.AddLine('x',[start[0],stop[0]])
        mesh.AddLine('y',[start[1],stop[1]])

        mesh.SmoothMeshLines('x',max_res/2)
        mesh.SmoothMeshLines('y',max_res/2)


        # apply the excitation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        start=[self.x-0.25, self.y-0.25, self.h-self.hs]
        stop =[self.x+0.25, self.y+0.25, self.h+self.hs]
        port= FDTD.AddLumpedPort(priority=5,port_nr=self.id,R=150,start=start,stop=stop,p_dir='z',excite=excite)
        return [CSX, FDTD, mesh, port]

