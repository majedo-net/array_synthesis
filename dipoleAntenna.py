from antenna import Antenna
from CSXCAD import ContinuousStructure
from openEMS import openEMS
from openEMS.physical_constants import *

class DipoleAntenna(Antenna):
    def __init__(self,id_,L_,W_,x_,y_,hs_,epsr_):
        super().__init__(id=id_,ant_type='Dipole')