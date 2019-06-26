from __future__ import absolute_import
from builtins import object
from proteus import *
from proteus.default_p import *
try:
    from .cylinder import *
except:
    from cylinder import *
from proteus.mprans import RANS3PF
name="rans3p"
LevelModelType = RANS3PF.LevelModel
if useOnlyVF:
    LS_model = None
else:
    LS_model = 2
if useRANS >= 1:
    Closure_0_model = 5; Closure_1_model=6
    if useOnlyVF:
        Closure_0_model=2; Closure_1_model=3
    if movingDomain:
        Closure_0_model += 1; Closure_1_model += 1
else:
    Closure_0_model = None
    Closure_1_model = None

coefficients = RANS3PF.Coefficients(epsFact=epsFact_viscosity,
                                    sigma=0.0,
                                    rho_0 = rho_0,
                                    nu_0 = nu_0,
                                    rho_1 = rho_1,
                                    nu_1 = nu_1,
                                    g=g,
                                    nd=nd,
                                    ME_model=V_model,
                                    PRESSURE_model=PRESSURE_model,
                                    SED_model=SED_model,
                                    VOF_model=VOF_model,
                                    VOS_model=VOS_model,
                                    LS_model=LS_model,
                                    Closure_0_model=Closure_0_model,
                                    Closure_1_model=Closure_1_model,
                                    epsFact_density=epsFact_density,
                                    stokes=False,
                                    useVF=useVF,
                                    useRBLES=useRBLES,
                                    useMetrics=useMetrics,
                                    eb_adjoint_sigma=1.0,
                                    eb_penalty_constant=weak_bc_penalty_constant,
                                    forceStrongDirichlet=ns_forceStrongDirichlet,
                                    turbulenceClosureModel=ns_closure,
                                    movingDomain=movingDomain,
                                    dragAlpha=dragAlpha,
                                    USE_SUPG=0.0,
                                    PSTAB=0.0,
                                    nParticles=1,
                                    particle_epsFact=1.5,
                                    particle_alpha=1000.0,
                                    particle_beta=1000.0,
                                    particle_penalty_constant=10.0,
                                    particle_sdfList=[particle_sdf],
                                    particle_velocityList=[particle_vel],
                                    use_ball_as_particle=1,
                                    ball_center=numpy.array([[0.2,0.2,0.0]]),
                                    ball_radius=numpy.array([0.05]),
                                    ball_velocity=numpy.array([[0.0,0.0,0.0]]),
                                    ball_angular_velocity=numpy.array([[0.0,0.0,0.0]]))
def zero(x, t):
    return 0.0
eps=1.0e-8
def getPeriodicBC(x,flag):
    if (x[0] < -L[0]/2.0+eps or x[0] >= L[0]/2.0-eps):
        return numpy.array([0.0,
                            round(x[1],5),
                            0.0])
def getDBC_u(x,flag):
    if flag in [boundaryTags['top'],boundaryTags['bottom']]:
        return lambda x,t: 0.0

def getDBC_v(x,flag):
    if flag in [boundaryTags['top'],boundaryTags['bottom']]:
        return lambda x,t: 0.0

def getAFBC_u(x,flag):
    if flag in [boundaryTags['left'],boundaryTags['right']]:
        return lambda x,t: 0.0

def getAFBC_v(x,flag):
    if flag in [boundaryTags['left'],boundaryTags['right']]:
        return lambda x,t: 0.0

def getDFBC_u(x,flag):
    if flag in [boundaryTags['left'],boundaryTags['right']]:
        return lambda x,t: 0.0

def getDFBC_v(x,flag):
    if flag in [boundaryTags['left'],boundaryTags['right']]:
        return lambda x,t: 0.0


advectiveFluxBoundaryConditions =  {0:getAFBC_u,
                                    1:getAFBC_v}

diffusiveFluxBoundaryConditions = {0:{0:getDFBC_u},
                                   1:{1:getDFBC_v}}


periodicDirichletConditions = {
                               0:getPeriodicBC,
                               1:getPeriodicBC}

class AtRest(object):
    def __init__(self):
        pass
    def uOfXT(self,x,t):
        return 0.0

initialConditions = {0:AtRest(),
                     1:AtRest()}
