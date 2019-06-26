from proteus import *
from proteus.default_p import *
from cylinder import *
from proteus.mprans import RANS2P

name = "momentum"

bcsTimeDependent = True

LevelModelType = RANS2P.LevelModel

coefficients = RANS2P.Coefficients(epsFact=epsFact_viscosity,
                                     rho_0=rho_0,
                                     nu_0=nu_0,
                                     rho_1=rho_1,
                                     nu_1=nu_1,
                                     g=g,
                                     nd=nd,
                                     LS_model=None,
                                     epsFact_density=epsFact_density,
                                     stokes=False,#useStokes,
                                     forceStrongDirichlet=ct.forceStrongDirichlet,
                                     eb_adjoint_sigma=1.0,
                                     eb_penalty_constant=10.0,
                                     useRBLES=0.0,
                                     useMetrics=1.0,
                                     use_ball_as_particle=use_ball_as_particle,
                                     ball_center=ball_center,
                                     ball_radius=ball_radius,
                                     ball_velocity=ball_velocity,
                                     ball_angular_velocity=ball_angular_velocity,
                                     nParticles = nParticles,
                                     particle_epsFact=1.5,
                                     particle_alpha=1e9,
                                     particle_beta=1e9,
                                     particle_penalty_constant=1e6,
                                     NONCONSERVATIVE_FORM = ct.nonconservative,
                                     MOMENTUM_SGE=ct.use_supg,
                                     PRESSURE_SGE=ct.use_supg,
                                     VELOCITY_SGE=ct.use_supg,
                                     PRESSURE_PROJECTION_STABILIZATION=0.0,
                                     nullSpace="NavierStokesConstantPressure")
#===============================================================================
# BC
#===============================================================================
def zero(x, t):
    return 0.0
eps=1.0e-8
def getPeriodicBC(x,tag):
    if (x[0] < -L[0]/2.0+eps or x[0] >= L[0]/2.0-eps):
        return numpy.array([0.0,
                            round(x[1],5),
                            0.0])

def getDBC_p(x,flag):
    return None

def getDBC_u(x,flag):
    if flag in [boundaryTags['top'],boundaryTags['bottom']]:
        return lambda x,t: 0.0

def getDBC_v(x,flag):
    if flag in [boundaryTags['top'],boundaryTags['bottom']]:
        return lambda x,t: 0.0

def getAFBC_p(x,flag):
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

dirichletConditions = {0:getDBC_p,
                       1:getDBC_u,
                       2:getDBC_v}

advectiveFluxBoundaryConditions =  {0:getAFBC_p,
                                    1:getAFBC_u,
                                    2:getAFBC_v}

diffusiveFluxBoundaryConditions = {0:{},
                                   1:{1:getDFBC_u},
                                   2:{2:getDFBC_v}}

periodicDirichletConditions = {0:getPeriodicBC,
                               1:getPeriodicBC,
                               2:getPeriodicBC}

class init_p:
    def uOfXT(self,x,t):
        return 0.0

class init_u:
    def uOfXT(self,x,t):
        return 0.0

class init_v:
    def uOfXT(self,x,t):
        return 0.0



initialConditions = {0:init_p(),
                     1:init_u(),
                     2:init_v()}
