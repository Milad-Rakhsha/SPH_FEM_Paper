from math import *
import proteus.MeshTools
from proteus import *
from proteus import Domain
from proteus.default_n import *
from proteus.Profiling import logEvent



L = (3.0, 1.5)
plate_dim=(0.02,0.3,0.02) # "Dimensions of the plate (Height/Width/thickness)"
plate_cent=(1.95,0.15,0.0) #Position of the center of the plate"),
plate_prop=(8000.0,5e6,0.3) #Physical Properties of the flexible plate (rho/E/nu)"),
plate_mesh_div=(1,30,2) #number of elements in each direction"),
dT_Chrono=0.0005

# Gravity
g = [0.0, -9.8]

# Initial condition
waterLine_x = 1
waterLine_z = 1.0

#  Discretization -- input options
#Refinement = 20#45min on a single core for spaceOrder=1, useHex=False
Refinement = 2
# Domain and mesh
he = L[0]/15.0
he*=0.25

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

# Time stepping/
T=8.0
dt_fixed = 0.01#0.03
dt_init = 0.0005 #min(0.1*dt_fixed,0.001)
runCFL=0.05
nDTout = int(round(T/dt_fixed))
tnList = [0.0,dt_init]+[dt_init+ i*dt_fixed for i in range(1,nDTout+1)]
info = open("TimeList.txt","w")


# Water
rho_0 = 1000.0
nu_0 = 1.0e-6

# Air
rho_1 = 1.205
nu_1 = 1.500e-5

# Sediment

rho_s = plate_prop[0]
nu_s = 1000.0*nu_0
dragAlpha = 0.0

# Surface tension
sigma_01 = 0.0



wall_x = 1.5
Um = 0.2



def signedDistance(x):
    phi_x = x[0] - waterLine_x
    phi_z = x[1] - waterLine_z
    if phi_x < 0.0:
        if phi_z < 0.0:
            return max(phi_x, phi_z)
        else:
            return phi_z
    else:
        if phi_z < 0.0:
            return phi_x
        else:
            return sqrt(phi_x ** 2 + phi_z ** 2)




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
        refined_L=max(plate_dim[0],plate_dim[1],plate_dim[2])*1.0
        vertices = [[0.0, 0.0],  #0
                    [L[0], 0.0], #1
                    [L[0], L[1]],#2
                    [0.0, L[1]], #3
                    # the following are set for refining the mesh
                    [plate_cent[0]-refined_L,0.0],
                    [plate_cent[0]+refined_L,0.0],
                    [plate_cent[0]+refined_L,plate_cent[1]+refined_L*0.8],
                    [plate_cent[0]-refined_L,plate_cent[1]+refined_L*0.8]
                    ]
        vertexFlags = [boundaryTags['left'],
                       boundaryTags['right'],
                       boundaryTags['right'],
                       boundaryTags['left'],
                        # the interior vertices should be flaged to 0
                       0, 0, 0, 0
                       ]

        segments = [[0, 1],
                    [1, 2],
                    [2, 3],
                    [3, 0],
                    #Interior segments
                    [4, 5],
                    [5, 6],
                    [6, 7],
                    [7, 4]
                    ]
        segmentFlags = [boundaryTags['bottom'],
                        boundaryTags['right'],
                        boundaryTags['top'],
                        boundaryTags['left'],
                        0,
                        0,
                        0,
                        0,
                        ]
        regions = [ [0.01,0.01],
                    [plate_cent[0],plate_cent[1]]
                   ]
        regionFlags = [1,
                       2
                       ]
        regionConstraints=[
            0.5*he**2,
            0.5*(he/2.0)**2
         ]


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

import Chrono_Proteus_Flex

class FlexiblePlate(AuxiliaryVariables.AV_base):
    def __init__(self,is3D=True,timeStep=dT_Chrono,
                m_plate_center=plate_cent,
                m_plate_dims=plate_dim,
                m_plate_num_div=plate_mesh_div,
                m_plate_prop=plate_prop,
                m_gravity=(0.0,0,0),
                m_free_x=(0,0,0),
                he=1.0,cfl_target=0.33,dt_init=dt_init):

        self.dt_init=dt_init
        self.chplate = Chrono_Proteus_Flex.FlexPlate(
                                        m_is3D=is3D,
                                        m_timeStep=timeStep,
                                        m_plate_center=np.array(m_plate_center,dtype="d"),
                                        m_plate_dims=np.array(m_plate_dims,dtype="d"),
                                        m_plate_num_div=np.array(m_plate_num_div,dtype="i"),
                                        m_plate_prop=np.array(m_plate_prop,dtype="d"),
                                        m_gravity=np.array(m_gravity,dtype="d"),
                                        m_free_x = np.array(m_free_x,dtype="i")
                                        )
        tmp_nnodes = self.chplate.getNumNodes()
        self.solidNodes = np.zeros((tmp_nnodes,3), 'd')
        self.solidNormals = np.zeros((tmp_nnodes,3), 'd')
        self.solidForces = np.zeros((tmp_nnodes, 3), 'd')
        self.new_coord=np.zeros((tmp_nnodes,3), 'd')
        self.new_coord_vel=np.zeros((tmp_nnodes,3), 'd')

    def attachModel(self,model,ar):
        self.chplate.attachModel(model,ar)
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


    def get_u(self):
        return 0.
    def get_v(self):
        return 0.
    def get_w(self):
        return 0.
    def calculate_init(self):
        self.last_F = None
        self.calculate()
    def getLocalNearestNode(self, location):
        # determine local nearest node distance
        nearest_node_distance_kdtree, nearest_node_kdtree = self.fluidNodes_kdtree.query(location)
        comm = Comm.get().comm.tompi4py()
        return comm.rank, nearest_node_kdtree, nearest_node_distance_kdtree

    def getLocalElement(self, femSpace, location, node):
        """Given a location and its nearest node, determine if it is on a
        local element.

        Returns None if location is not on any elements owned by this
        process

        """

        # search elements that contain the nearest node
        patchBoundaryNodes=set()
        checkedElements=[]
        for eOffset in range(femSpace.mesh.nodeElementOffsets[node], femSpace.mesh.nodeElementOffsets[node + 1]):
            eN = femSpace.mesh.nodeElementsArray[eOffset]
            checkedElements.append(eN)
            patchBoundaryNodes|=set(femSpace.mesh.elementNodesArray[eN])
            # evaluate the inverse map for element eN
            xi = femSpace.elementMaps.getInverseValue(eN, location)
            # query whether xi lies within the reference element
            if femSpace.elementMaps.referenceElement.onElement(xi):
                return eN
        for node in patchBoundaryNodes:
            for eOffset in range(femSpace.mesh.nodeElementOffsets[node], femSpace.mesh.nodeElementOffsets[node + 1]):
                eN = femSpace.mesh.nodeElementsArray[eOffset]
                if eN not in checkedElements:
                    checkedElements.append(eN)
                    # evaluate the inverse map for element eN
                    xi = femSpace.elementMaps.getInverseValue(eN, location)
                    # query whether xi lies within the reference element
                    if femSpace.elementMaps.referenceElement.onElement(xi):
                        return eN
        # no elements found
        return None

    def findNearestNode(self, femSpace, location):
        """Given a gauge location, attempts to locate the most suitable
        process for monitoring information about this location, as
        well as the node on the process closest to the location.

        Returns a 2-tuple containing an identifier for the closest
        'owning' process as well as the local ids of the node and
        nearest element.

        """
        from mpi4py import MPI
        comm = Comm.get().comm.tompi4py()
        comm_rank, nearest_node, nearest_node_distance = self.getLocalNearestNode(location)
        local_element = self.getLocalElement(femSpace, location, nearest_node)

        # determine global nearest node
        haveElement = int(local_element is not None)
        global_have_element, owning_proc = comm.allreduce((haveElement, comm.rank),
                                                          op=MPI.MAXLOC)
        if global_have_element:
            # logEvent("Gauges on element at location: [%g %g %g] assigned to %d" % (location[0], location[1], location[2],
            #                                                         owning_proc), 3)
            pass
        else:
            # gauge isn't on any of the elements, just use nearest node
            global_min_distance, owning_proc = comm.allreduce((nearest_node_distance, comm.rank), op=MPI.MINLOC)
            # logEvent("Off-element gauge location: [%g %g %g] assigned to %d" % (location[0], location[1], location[2],
            #                                                      owning_proc), 3)
        if comm.rank != owning_proc:
            nearest_node = None
        # logEvent("owning_proc "+`owning_proc`)
        # logEvent("nearest_node "+`nearest_node`)
        assert owning_proc is not None
        return owning_proc, nearest_node

    def calculate(self):
        import  numpy as np
        from numpy.linalg import inv
        import copy
        self.solidNodes, self.solidNormals = self.chplate.calcNodalInfo()

        size=self.solidNodes.shape[0]
        #print(self.solidNodes)
        # print(self.solidNormals)

        from mpi4py import MPI
        from scipy import spatial
        #initialy put the whole point finding and calculation algorithm here
        #recompute kdtree if mesh is moving
        comm = Comm.get().comm.tompi4py()
        rank = comm.Get_rank()
        self.fluidNodes_kdtree = spatial.cKDTree(self.model.levelModelList[-1].mesh.nodeArray)
        # self.solidForces[:,:] = 0.0
        self.solidForces=np.zeros((size,3))
        # if (rank==0):
        #     print self.solidNodes, self.solidNormals
        for i, x, n in zip(range(self.solidNodes.shape[0]), self.solidNodes, self.solidNormals):
            rank_of_node, node = self.findNearestNode(self.model.levelModelList[-1].u[0].femSpace, x)
            if rank_of_node == comm.rank:
                eN = self.getLocalElement(self.model.levelModelList[-1].u[0].femSpace, x, node)
                if eN is None:
                    p = self.model.levelModelList[-1].pressureModel.u[0].dof[node]
                    u = self.model.levelModelList[-1].u[0].dof[node]
                    v = self.model.levelModelList[-1].u[1].dof[node]
                    #gradients for viscous stress
                    if self.nd > 2:
                        w = self.model.levelModelList[-1].u[3].dof[node]
                    if self.nd <= 2:
                        w = 0
                        # grad_u = self.model.levelModelList[-1].u[1].getGradientValue(eN, xi)
                else:
                    # logEvent("i, x, n, eN "+`i`+","+`x`+","+`n`+","+`eN`)
                    xi = self.model.levelModelList[-1].u[0].femSpace.elementMaps.getInverseValue(eN, x)
                    p = self.model.levelModelList[-1].pressureModel.u[0].getValue(eN,xi)
                    u = self.model.levelModelList[-1].u[0].getValue(eN, xi)
                    v = self.model.levelModelList[-1].u[1].getValue(eN, xi)
                    #gradients for viscous stress
                    if self.nd > 2:
                        w = self.model.levelModelList[-1].u[3].getValue(eN, xi)
                    if self.nd <= 2:
                        w = 0

                # grad_p= self.model.levelModelList[-1].pressureModel.u[0].getGradientValue(eN, xi)
                # grad_u = self.model.levelModelList[-1].u[0].getGradientValue(eN, xi)
                # grad_v = self.model.levelModelList[-1].u[1].getGradientValue(eN, xi)
                #if self.nd > 2:
                #    grad_w = self.model.levelModelList[-1].u[3].getGradientValue(eN, xi)
                #if self.nd <= 2:
                #    grad_w = np.zeros((3,),'d')
                 #print(self.model.levelModelList )
                #vof = self.model.levelModelList[1].u[0].dof[node]

                #H_rho = useVF*min(1.0,max(0.0,vof))
                #rho = rho_0*(1.0-H_rho)+rho_1*H_rho
                # print(p,n)
                self.solidForces[i, 0] = -p*n[0] #+ mu*(grad_u[0][0] + grad_u[0][0])*n[0] + mu*(grad_u[0][1] + grad_u[1][0])*n[1] + mu*(grad_u[0][2] + grad_u[2][0])*n[2]
                self.solidForces[i, 1] = -p*n[1] #+ mu*(grad_u[0][1] + grad_u[1][0])*n[0] + mu*(grad_u[1][1] + grad_u[1][1])*n[1] + mu*(grad_u[1][2] + grad_u[2][1])*n[2]
                self.solidForces[i, 2] = -p*n[2] #+ mu*(grad_u[0][2] + grad_u[2][0])*n[0] + mu*(grad_u[2][1] + grad_u[1][2])*n[1] + mu*(grad_u[2][2] + grad_u[2][2])*n[2]
        #cek hack, can do this communication more efficiently
        comm.Allreduce([self.solidForces.copy(), MPI.DOUBLE],
                       [self.solidForces, MPI.DOUBLE],MPI.SUM)

        comm = MPI.COMM_WORLD
        nprocs = comm.Get_size()
        rank   = comm.Get_rank()
        try:
            self.proteus_dt = self.model.levelModelList[-1].dt_last
            t = self.model.stepController.t_model_last
        except:
            self.proteus_dt = self.dt_init
            t = 0


        if (t>=self.dt_init):
            if rank==0:
                print("time/dt before chrono calculation:" , t , self.proteus_dt)
                # print self.solidForces
                self.new_coord,self.new_coord_vel=self.chplate.calculate(self.solidForces,self.proteus_dt)
                for writeTime in tnList:
                    if (t>0.0001):
                        if (abs(t-writeTime)<1e-6):
                            self.chplate.writeFrame()

            self.new_coord = comm.bcast(self.new_coord, root=0)
            self.new_coord_vel=comm.bcast(self.new_coord_vel, root=0)
            comm.Barrier()

            # SyncData processes the new data by calling "prepareData" on the chrono object
            if rank!=0:
                # print ("processor ", rank)
                self.chplate.SyncData(self.new_coord,self.new_coord_vel)
                # self.chplate.calcNodalInfo()
        else:
                self.chplate.calcNodalInfo()

plate = FlexiblePlate(is3D=False,timeStep=dT_Chrono,m_plate_center=plate_cent,
                      m_plate_dims=plate_dim,m_plate_num_div=plate_mesh_div,
                      m_plate_prop=plate_prop,m_gravity=(0.0,0.0,0.0),m_free_x=(0,0,0),
                      he=1.0,cfl_target=0.1,dt_init=dt_init)

def particle_sdf(t, x):
    N=np.zeros((3,1), 'd')
    d , N=plate.chplate.d_N_IBM(x)
    return d-0.0,(N[0],N[1])

import numpy as np
def particle_vel(t, x):
    v=np.zeros((3,1), 'd')
    v=plate.chplate.vel_IBM(x)

    return (v[0], v[1])
