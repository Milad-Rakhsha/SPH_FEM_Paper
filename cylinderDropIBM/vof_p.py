from proteus import *
from proteus.default_p import *
from proteus.ctransportCoefficients import smoothedHeaviside
from cylinderDrop import *
from proteus.mprans import VOF3P

LevelModelType = VOF3P.LevelModel

coefficients = VOF3P.Coefficients(LS_model=LS_model,
                                  V_model=V_model,
                                  RD_model=RD_model,
                                  ME_model=VOF_model,
                                  VOS_model=VOS_model,
                                  checkMass=True,
                                  useMetrics=useMetrics,
                                  epsFact=epsFact_vof,
                                  sc_uref=vof_sc_uref,
                                  sc_beta=vof_sc_beta,
                                  movingDomain=movingDomain)

def getDBC_vof(x,flag):
    if flag in [boundaryTags['left'],boundaryTags['right']]:
        return lambda x,t: smoothedHeaviside(epsFact_consrv_heaviside*he,signedDistance(x))

dirichletConditions = {0:getDBC_vof}

def getAFBC_vof(x,flag):
    if flag not in [boundaryTags['left'],boundaryTags['right']]:
        return lambda x,t: 0.0

advectiveFluxBoundaryConditions = {0:getAFBC_vof}
diffusiveFluxBoundaryConditions = {0:{}}

class VOF_IC:
    def uOfXT(self,x,t):
        return smoothedHeaviside(epsFact_consrv_heaviside*he,signedDistance(x))

initialConditions  = {0:VOF_IC()}
