
import numpy as np
from array_funcs import *
from circ_rps import *
from entry import *
import matplotlib.pyplot as plt
from matplotlib.patches import Circle

plt.rc('font', family='serif')

if __name__ == '__main__':

    # Optimized Array
    h= 33 

    xs= [ 0.00000000e+00,  6.52103245e-02,  3.34365927e-02, -3.17737318e-02,
 -6.52103245e-02, -3.34365927e-02,  3.17737318e-02,  1.82159034e-01,
  1.42392805e-01,  1.92148025e-02, -1.15218971e-01, -1.82159034e-01,
 -1.42392805e-01, -1.92148025e-02,  1.15218971e-01,  3.21942168e-01,
  2.05757139e-01,  1.09798765e-02, -1.87991325e-01, -3.15156231e-01,
 -3.21942168e-01, -2.05757139e-01, -1.09798765e-02,  1.87991325e-01,
  3.15156231e-01,  5.14450540e-01,  4.45527237e-01,  2.57225270e-01,
  3.15010104e-17, -2.57225270e-01, -4.45527237e-01, -5.14450540e-01,
 -4.45527237e-01, -2.57225270e-01, -9.45030311e-17,  2.57225270e-01,
  4.45527237e-01] 

    ys= [ 0.00000000e+00, -9.60053182e-04,  5.59937710e-02,  5.69538242e-02,
  9.60053182e-04, -5.59937710e-02, -5.69538242e-02, -1.92148025e-02,
  1.15218971e-01,  1.82159034e-01,  1.42392805e-01,  1.92148025e-02,
 -1.15218971e-01, -1.82159034e-01, -1.42392805e-01,  9.30604259e-02,
  2.64520324e-01,  3.34942449e-01,  2.77427943e-01,  1.13945392e-01,
 -9.30604259e-02, -2.64520324e-01, -3.34942449e-01, -2.77427943e-01,
 -1.13945392e-01,  0.00000000e+00,  2.57225270e-01,  4.45527237e-01,
  5.14450540e-01,  4.45527237e-01,  2.57225270e-01,  6.30020208e-17,
 -2.57225270e-01, -4.45527237e-01, -5.14450540e-01, -4.45527237e-01,
 -2.57225270e-01] 

    spirads= [0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.06, 0.06, 0.06, 0.06, 0.06, 0.06, 0.06,
 0.06, 0.09, 0.08, 0.07, 0.07, 0.08, 0.09, 0.08, 0.07, 0.07, 0.08, 0.1,  0.1,  0.09,
 0.09, 0.09, 0.09, 0.1,  0.1,  0.09, 0.09, 0.09, 0.09]


    fmax = 2.3e9
    fmin = 0.4e9
    freq = fmax * (1+np.sin(60*np.pi/180))
    element_pattern_max = fetch_element_patterns(spirads,freq,h)
    element_pattern_min = fetch_element_patterns(spirads,fmin,h)
    lam_max = 3e8/fmax
    lam_min = 3e8/fmin
    k_max = 2*np.pi/lam_max
    k_min = 2*np.pi/lam_min
    theta = np.linspace(-np.pi/2, np.pi/2, 181)
    phi = np.linspace(0, np.pi, 181)

    # make scanned plots
    theta_0 = np.deg2rad(60)
    phi_0 = np.deg2rad(135)

    ArrF_max, Tot_max = scan_array_factor(xs,ys,k_max,element_pattern_max,theta,phi,theta_0,phi_0)
    ArrF_min, Tot_min = scan_array_factor(xs,ys,k_min,element_pattern_min,theta,phi,theta_0,phi_0)
    ArrF_max = (25/35)*ArrF_max

    makePatternPlotsOnlyTheta(theta,phi,ArrF_max,Tot_max,0,fmax,'optimized',save=True)
    makePatternPlotsOnlyTheta(theta,phi,ArrF_min,Tot_min,0,fmin,'optimized',save=True)
