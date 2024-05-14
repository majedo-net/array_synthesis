from array_funcs import makeUVPlot
from os import listdir
from os.path import join
import numpy as np

filelist = [f for f in listdir('results/aws-tess-array/nr10') if 'ff' in f]

theta = np.linspace(0, np.pi, 181)
phi = np.linspace(0, 2*np.pi, 361)

for f in filelist:
    fname = f.rstrip('.txt')
    fname = 'results/aws-tess-array/nr10/'+fname+'.png'
    G = np.genfromtxt(join('results/aws-tess-array/nr10',f))
    G = 10*np.log10(G)
    makeUVPlot(theta,phi,G,'Embedded Pattern',fname)

print('breakpoint')
