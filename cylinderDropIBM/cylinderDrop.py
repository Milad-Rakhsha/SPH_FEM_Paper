from math import *
import proteus.MeshTools
from proteus import *
from proteus import Domain
from proteus.default_n import *
from proteus.Profiling import logEvent



g_chrono = [0.0, -9.81, 0.0]
g = [0.0, -9.81, 0.0]

# Domain and mesh
L = [1.0,1.6,0.55]
# Initial condition
container_dim=[L[0],L[1],L[2]] #Dimensions of the container (Height/Width/depth)
container_cent=[L[0]/2,L[1]/2,0.0] #Position of the center of the container"
#Dimensions of the particles in initial configuration
# (Height/Width/thickness) in case you have a box of particles
# (Radius,unused, unused) if you have a spherical set of particles
particle_diameter=0.24
particle_density=1000.0*0.7
particle_dim=[particle_diameter,0.2,0.0] #Radius/Height
particle_cent=[L[0]/2, 1.25, 0.0] #center of the container (Height/Width/depth)
dT_Chrono=0.001



# Time stepping/
T=10.0
dt_fixed = 0.01#0.03
dt_init = 0.001 #min(0.1*dt_fixed,0.001)
runCFL=0.05
nDTout = int(round(T/dt_fixed))
tnList = [0.0,dt_init]+[i*dt_fixed for i in range(1,nDTout+1)]




nLevels = 1
#parallelParti

# Water
rho_0 = 1000.0
nu_0 = 1.0e-6

# Air
rho_1 = 1.205
nu_1 = 1.500e-5

# Sediment

rho_s = rho_0
nu_s = 10000.0*nu_0
dragAlpha = 0.0

# Surface tension
sigma_01 = 0.0

# Initial condition
waterLine_x = 0.75
waterLine_z = 0.2



#  Discretization -- input options
#Refinement = 20#45min on a single core for spaceOrder=1, useHex=False
#  Discretization -- input options
#Refinement = 20#45min on a single core for spaceOrder=1, useHex=False
Refinement = 2
# Domain and mesh
he = 0.025

sedimentDynamics=False
genMesh = True
movingDomain = False
applyRedistancing = True
useOldPETSc = False
useSuperlu = False
timeDiscretization = 'vbdf'#vbdf'#'vbdf'  # 'vbdf', 'be', 'flcbdf'
spaceOrder = 2
pspaceOrder = 1
useHex = False
useRBLES = 0.0
useMetrics = 1.0
applyCorrection = True
useVF = 0.0
useOnlyVF = False
useRANS = 0  #  -- None
             # 1 -- K-Epsilon
             # 2 -- K-Omega
openTop=True
fl_H = L[1]



nLevels = 1
#parallelPartitioningType = proteus.MeshTools.MeshParallelPartitioningTypes.element
parallelPartitioningType = proteus.MeshTools.MeshParallelPartitioningTypes.node
nLayersOfOverlapForParallel = 0
structured = False


# Input checks
if spaceOrder not in [1, 2]:
    print "INVALID: spaceOrder" + spaceOrder
    sys.exit()

if useRBLES not in [0.0, 1.0]:
    print "INVALID: useRBLES" + useRBLES
    sys.exit()

if useMetrics not in [0.0, 1.0]:
    print "INVALID: useMetrics"
    sys.exit()

#  Discretization
nd = 2

if spaceOrder == 1:
    hFactor = 1.0
    if useHex:
        basis = C0_AffineLinearOnCubeWithNodalBasis
        elementQuadrature = CubeGaussQuadrature(nd, 2)
        elementBoundaryQuadrature = CubeGaussQuadrature(nd - 1, 2)
    else:
        basis = C0_AffineLinearOnSimplexWithNodalBasis
        elementQuadrature = SimplexGaussQuadrature(nd, 3)
        elementBoundaryQuadrature = SimplexGaussQuadrature(nd - 1, 3)
elif spaceOrder == 2:
    hFactor = 0.5
    if useHex:
        basis = C0_AffineLagrangeOnCubeWithNodalBasis
        elementQuadrature = CubeGaussQuadrature(nd, 4)
        elementBoundaryQuadrature = CubeGaussQuadrature(nd - 1, 4)
    else:
        basis = C0_AffineQuadraticOnSimplexWithNodalBasis
        elementQuadrature = SimplexGaussQuadrature(nd, 5)
        elementBoundaryQuadrature = SimplexGaussQuadrature(nd - 1, 5)

if pspaceOrder == 1:
    if useHex:
        pbasis = C0_AffineLinearOnCubeWithNodalBasis
    else:
        pbasis = C0_AffineLinearOnSimplexWithNodalBasis
elif pspaceOrder == 2:
    if useHex:
        pbasis = C0_AffineLagrangeOnCubeWithNodalBasis
    else:
        pbasis = C0_AffineQuadraticOnSimplexWithNodalBasis





# new style PointGauges
# pointGauges = PointGauges(gauges = ((('u', 'v'), ((0.5, 0.5, 0), (1, 0.5, 0))), (('p',), ((0.5, 0.5, 0),))),
#                           activeTime=(0, 0.5),
#                           sampleRate=0,
#                           fileName='combined_gauge_0_0.5_sample_all.csv')

# lineGauges = LineGauges(gaugeEndpoints={'lineGauge_xtoH=0.825': ((0.495, 0.0, 0.0), (0.495, 1.8, 0.0))}, linePoints=20)
# #'lineGauge_x/H=1.653':((0.99,0.0,0.0),(0.99,1.8,0.0))
# lineGauges_phi = LineGauges_phi(lineGauges.endpoints, linePoints=20)

if useHex:
    nnx = 4 * Refinement + 1
    nny = 2 * Refinement + 1
    hex = True
    domain = Domain.RectangularDomain(L)
else:
    boundaries = ['left', 'right', 'bottom', 'top', 'front', 'back']
    boundaryTags = dict([(key, i + 1) for (i, key) in enumerate(boundaries)])
    if structured:
        nnx = 4 * Refinement
        nny = 2 * Refinement
    else:
        vertices = [[0.0, 0.0],  #0
                    [L[0], 0.0], #1
                    [L[0], L[1]],#2
                    [0.0, L[1]]] #3]
        vertexFlags = [boundaryTags['left'],
                       boundaryTags['right'],
                       boundaryTags['right'],
                       boundaryTags['left']]
                        # the interior vertices should be flaged to 0


        segments = [[0, 1],
                    [1, 2],
                    [2, 3],
                    [3, 0]
                    ]
        segmentFlags = [boundaryTags['bottom'],
                        boundaryTags['right'],
                        boundaryTags['top'],
                        boundaryTags['left']]

        regions = [[0.01,0.01] ]
        regionFlags = [1]
        regionConstraints=[0.5*he**2]


        #        for gaugeName,gaugeCoordinates in pointGauges.locations.iteritems():
        #            vertices.append(gaugeCoordinates)
        #            vertexFlags.append(pointGauges.flags[gaugeName])

        # for gaugeName, gaugeLines in lineGauges.linepoints.iteritems():
        #     for gaugeCoordinates in gaugeLines:
        #         vertices.append(gaugeCoordinates)
        #         vertexFlags.append(lineGauges.flags[gaugeName])
        domain = Domain.PlanarStraightLineGraphDomain(vertices=vertices,
                                                      vertexFlags=vertexFlags,
                                                      segments=segments,
                                                      segmentFlags=segmentFlags,
                                                      regions=regions,
                                                      regionFlags=regionFlags,
                                                      regionConstraints=regionConstraints)
        #go ahead and add a boundary tags member
        domain.boundaryTags = boundaryTags
        domain.writePoly("mesh")
        domain.writePLY("mesh")
        domain.writeAsymptote("mesh")
        # triangleOptions = "VApq30Dena%8.8f" % ((he ** 2) / 2.0,)
        triangleOptions = "VApq30Dena"

logEvent("""Mesh generated using: tetgen -%s %s""" % (triangleOptions, domain.polyfile + ".poly"))


# Numerical parameters
ns_forceStrongDirichlet = False
ns_sed_forceStrongDirichlet = False
if useMetrics:
    ns_shockCapturingFactor  = 0.5
    ns_lag_shockCapturing = True
    ns_lag_subgridError = True
    ls_shockCapturingFactor  = 0.5
    ls_lag_shockCapturing = True
    ls_sc_uref  = 1.0
    ls_sc_beta  = 1.5
    vof_shockCapturingFactor = 0.5
    vof_lag_shockCapturing = True
    vof_sc_uref = 1.0
    vof_sc_beta = 1.5
    rd_shockCapturingFactor  = 0.5
    rd_lag_shockCapturing = False
    epsFact_density    = 1.5
    epsFact_viscosity  = epsFact_curvature  = epsFact_vof = epsFact_consrv_heaviside = epsFact_consrv_dirac = epsFact_density
    epsFact_redistance = 0.33
    epsFact_consrv_diffusion = 10.0
    redist_Newton = True
    kappa_shockCapturingFactor = 0.25
    kappa_lag_shockCapturing = True#False
    kappa_sc_uref = 1.0
    kappa_sc_beta = 1.0
    dissipation_shockCapturingFactor = 0.25
    dissipation_lag_shockCapturing = True#False
    dissipation_sc_uref = 1.0
    dissipation_sc_beta = 1.0
else:
    ns_shockCapturingFactor = 0.9
    ns_lag_shockCapturing = True
    ns_lag_subgridError = True
    ns_sed_shockCapturingFactor = 0.9
    ns_sed_lag_shockCapturing = True
    ns_sed_lag_subgridError = True
    ls_shockCapturingFactor = 0.9
    ls_lag_shockCapturing = True
    ls_sc_uref = 1.0
    ls_sc_beta = 1.0
    vof_shockCapturingFactor = 0.9
    vof_lag_shockCapturing = True
    vof_sc_uref = 1.0
    vof_sc_beta = 1.0
    vos_shockCapturingFactor = 0.9
    vos_lag_shockCapturing = True
    vos_sc_uref = 1.0
    vos_sc_beta = 1.0
    rd_shockCapturingFactor = 0.9
    rd_lag_shockCapturing = False
    epsFact_density = 1.5
    epsFact_viscosity = epsFact_curvature = epsFact_vof = epsFact_vos = epsFact_consrv_heaviside = epsFact_consrv_dirac = \
        epsFact_density
    epsFact_redistance = 0.33
    epsFact_consrv_diffusion = 0.1
    redist_Newton = False
    kappa_shockCapturingFactor = 0.9
    kappa_lag_shockCapturing = True  #False
    kappa_sc_uref = 1.0
    kappa_sc_beta = 1.0
    dissipation_shockCapturingFactor = 0.9
    dissipation_lag_shockCapturing = True  #False
    dissipation_sc_uref = 1.0
    dissipation_sc_beta = 1.0

ns_nl_atol_res = max(1.0e-10, 0.01 * he ** 2)
ns_sed_nl_atol_res = max(1.0e-10, 0.01 * he ** 2)
vof_nl_atol_res = max(1.0e-10, 0.01 * he ** 2)
vos_nl_atol_res = max(1.0e-10, 0.01 * he ** 2)
ls_nl_atol_res = max(1.0e-10, 0.01 * he ** 2)
rd_nl_atol_res = max(1.0e-10, 0.05 * he)
mcorr_nl_atol_res = max(1.0e-10, 0.01 * he ** 2)
kappa_nl_atol_res = max(1.0e-10, 0.01 * he ** 2)
dissipation_nl_atol_res = max(1.0e-10, 0.01 * he ** 2)
phi_nl_atol_res = max(1.0e-10, 0.01 * he ** 2)
pressure_nl_atol_res = max(1.0e-10, 0.01 * he ** 2)

#turbulence
ns_closure = 0  #1-classic smagorinsky, 2-dynamic smagorinsky, 3 -- k-epsilon, 4 -- k-omega
ns_sed_closure = 0  #1-classic smagorinsky, 2-dynamic smagorinsky, 3 -- k-epsilon, 4 -- k-omega
if useRANS == 1:
    ns_closure = 3
elif useRANS == 2:
    ns_closure == 4


print ('importing chrono model')
import Chrono
print('done!')
class ChronoModel(AuxiliaryVariables.AV_base):
    def __init__(self,timeStep=1e-3,
                m_container_center=container_cent,
                m_container_dims=container_dim,
                m_particles_center=particle_cent,
                m_particles_dims=particle_dim,
                m_gravity=g_chrono,
                m_particles_diameter=particle_diameter,
                m_particles_density=particle_density,
                dt_init=dt_init):
        self.mtime=0
        self.dt_init=dt_init
        self.chmodel = Chrono.MBDModel(
                                        m_timeStep=timeStep,
                                        m_container_center=np.array(m_container_center,dtype="d"),
                                        m_container_dims=np.array(m_container_dims,dtype="d"),
                                        m_particles_center=np.array(m_particles_center,dtype="d"),
                                        m_particles_dims=np.array(m_particles_dims,dtype="d"),
                                        m_particles_diameter=particle_diameter,
                                        m_particles_density=particle_density,
                                        m_gravity=np.array(m_gravity,dtype="d")
                                        )
        tmp_nnodes = self.chmodel.get_num_surface_points()
        self.solidNodes = np.zeros((tmp_nnodes,3), 'd')
        self.solidNormals = np.zeros((tmp_nnodes,3), 'd')
        self.solidForces = np.zeros((tmp_nnodes, 3), 'd')
        self.new_coord=np.zeros((tmp_nnodes,3), 'd')

    def attachModel(self,model,ar):
        self.chmodel.attachModel(model,ar)
        self.model=model
        self.ar=ar
        self.writer = Archiver.XdmfWriter()
        self.nd = model.levelModelList[-1].nSpace_global
        m = self.model.levelModelList[-1]
        flagMax = max(m.mesh.elementBoundaryMaterialTypes)
        flagMin = min(m.mesh.elementBoundaryMaterialTypes)
        assert(flagMin >= 0)
        assert(flagMax <= 7)
        self.nForces=flagMax+1
        assert(self.nForces <= 8)
        return self

    def calculate_init(self):
        self.last_F = None
        self.calculate()

    def calculate(self):
        import  numpy as np
        from numpy.linalg import inv
        import copy
        self.solidForces = self.model.levelModelList[-1].coefficients.particle_netForces;
        # self.solidMoments=self.model.levelModelList[-1].coefficients.particle_netForces;

        try:
            self.proteus_dt = self.model.levelModelList[-1].dt_last
            t = self.model.stepController.t_model_last
        except:
            self.proteus_dt = self.dt_init
            t = 0

        print("time/dt before chrono calculation:" , t , self.proteus_dt)
        print(self.solidForces)
        if (t>=self.dt_init):
            self.chmodel.calculate(self.solidForces, self.proteus_dt)
            print "done with chrono"
        else:
            f= open("Forces.txt","w+")
            f.write("t,Fx,Fy\n")
            f.close()
        for writeTime in tnList:
            if (t>0.0001):
                if (abs(t-writeTime)<0.0001):
                    self.chmodel.writeFrame()
                    f= open("Forces.txt","a+")
                    print self.solidForces[0]
                    f.write("%f,%f,%f\n" % (t,self.solidForces[0,0],self.solidForces[0,1]))
                    f.close()


myChModel = ChronoModel(timeStep=dT_Chrono,
                    m_container_center=container_cent,m_container_dims=container_dim,
                    m_particles_center=particle_cent,m_particles_dims=particle_dim,
                    m_particles_diameter=particle_diameter,m_particles_density=particle_density,
                    m_gravity=(g_chrono[0],g_chrono[1],g_chrono[2]))

def signedDistance(x):
    return x[1] - 1.0

def obj_sdf_Calc(t,x):
    d,n=myChModel.chmodel.calc_d_N_IBM(x)
    return d,(n[0],n[1])

def obj_vel_Calc(t,x):
    v=myChModel.chmodel.calc_v_IBM(x)
    return (v[0],v[1])
