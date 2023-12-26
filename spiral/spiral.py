import os, tempfile
from pylab import *
import numpy as np

from CSXCAD  import ContinuousStructure
from openEMS import openEMS
from openEMS.physical_constants import *

def makeSpiral():

    return
    

if __name__ == '__main__':
    centers = np.ndarray([[0, 0], [200, 100],[100, -100],[-100, -100],[-100, 100]])
    radii = np.ndarray([[5,50],[5,40],[5,30],[5,20]])

    freq = 2e9
    unit = 1e-3 # all length in mm
    f_start = 0.9 * freq
    f_stop = freq
    max_res = np.floor(C0 / (f_stop) / unit / 20) #cell size: lambda/30
    padding = max_res *20
    hs = 1.6 # substrate thickness
    h = 30 # cavity height
    phase_center = np.ndarray([np.mean(centers,axis=0), h])
