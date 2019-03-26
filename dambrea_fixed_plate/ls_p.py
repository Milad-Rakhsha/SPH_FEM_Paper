from proteus import *
from proteus.default_p import *
from dambreak import *
from proteus.mprans import NCLS

LevelModelType = NCLS.LevelModel

coefficients = NCLS.Coefficients(V_model=V_model,
                                   RD_model=RD_model,
                                   ME_model=LS_model,
                                   checkMass=False,
                                   useMetrics=useMetrics,
                                   epsFact=epsFact_consrv_heaviside,
                                   sc_uref=ls_sc_uref,
                                   sc_beta=ls_sc_beta,
                                   movingDomain=movingDomain)

def getDBC_ls(x,flag):
    return None

dirichletConditions = {0:getDBC_ls}

advectiveFluxBoundaryConditions =  {}
diffusiveFluxBoundaryConditions = {0:{}}

class PHI_IC:
    def uOfXT(self,x,t):
        return signedDistance(x)

initialConditions  = {0:PHI_IC()}
