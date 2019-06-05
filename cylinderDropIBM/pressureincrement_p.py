from math import *
from proteus import *
from proteus.default_p import *
from cylinderDrop import *


#domain = ctx.domain
#nd = ctx.nd
name = "pressureincrement"

#from ProjectionScheme import PressureIncrement
#coefficients=PressureIncrement(rho_f_min = rho_1,
#                               rho_s_min = rho_s,
#                               nd = nd,
#                               modelIndex=PINC_model,
#                               fluidModelIndex=V_model)
from proteus.mprans import PresInc
coefficients=PresInc.Coefficients(rho_f_min = (1.0-1.0e-8)*rho_1,
                                 rho_s_min = (1.0-1.0e-8)*rho_s,
                                 nd = nd,
                                 modelIndex=PINC_model,
                                 fluidModelIndex=V_model)

LevelModelType = PresInc.LevelModel

#pressure increment should be zero on any pressure dirichlet boundaries
def getDBC_phi(x,flag):
    if flag in [boundaryTags['top']]: #,boundaryTags['left'],boundaryTags['front'], boundaryTags['back']]:
        return lambda x,t: 0.0
#def getDBC_phi(x,flag):
#    if flag == boundaryTags['top'] and openTop:
#        return lambda x,t: 0.0
from math import cos,pi

#the advectiveFlux should be zero on any no-flow  boundaries
def getAdvectiveFlux_qt(x,flag):
    if flag in[ boundaryTags['left'],boundaryTags['right'],boundaryTags['front'],boundaryTags['back'],boundaryTags['bottom']]:
        return lambda x,t: 0.0
    else:
        return None

def getDiffusiveFlux_phi(x,flag):
    if flag == boundaryTags['top']:
        return None
    else:
        return lambda x,t: 0.0

class getIBC_phi:
    def __init__(self):
        pass
    def uOfXT(self,x,t):
        return 0.0

initialConditions = {0:getIBC_phi()}

dirichletConditions = {0:getDBC_phi }
advectiveFluxBoundaryConditions = {0:getAdvectiveFlux_qt}
diffusiveFluxBoundaryConditions = {0:{0:getDiffusiveFlux_phi}}
