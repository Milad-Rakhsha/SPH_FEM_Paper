import numpy
cimport numpy
from proteus import AuxiliaryVariables, Archiver
from proteus.Profiling import  logEvent
from cython.view cimport array as cvarray

import numpy as np
# Import the C-level symbols of numpy
cimport numpy as np
from libcpp cimport bool
cimport proteus.mbd.CouplingFSI as fsi
from libcpp.memory cimport (shared_ptr, make_shared)
# chrono C++ headers
cimport proteus.mbd.ChronoHeaders as ch

cdef extern from "Chrono_Proteus_Flex.h":
    cdef cppclass cppChFlexPlate:
        void step(double* forces, int num_nodes, double dt)
        double* calcNodalPos()
        double* calcNodalNormal()
        void calc_d_N_IBM(double* x, double* out)
        void calc_vel_IBM(double* x, double* out)
        void prepareData(double *nodal_corrd_v, double *nodal_corrd_vel_v)
        const int num_surface_nodes
        const int num_total_nodes
        double *nodal_corrd_vec
        double *nodal_corrd_vel_vec
        void writeThisFrame()


    cppChFlexPlate* newChFlexPlate(shared_ptr[ch.ChSystemSMC] _system,
                                    bool m_is3D,
                                    double m_timeStep,
                                    double* m_plate_center,
                                    double* m_plate_dims,
                                    int* m_plate_num_div,
                                    double* m_plate_prop,
                                    double* m_gravity,
                                    int* m_free_x)

cdef class FlexPlate:
    cdef cppChFlexPlate* thisptr
    cdef object model
    cdef object np_pos
    cdef object np_normal
    cdef object n_coor
    cdef object v_coor

    cdef object numnodes

    cdef fsi.ProtChSystem my_system

    def __cinit__(self,
                  fsi.ProtChSystem _system,
                  bool m_is3D=True,
                  double m_timeStep=5e-3,
                  numpy.ndarray m_plate_center=numpy.array((0.0,0.0,0.0),dtype="d"),
                  numpy.ndarray m_plate_dims=numpy.array((2.0,1.0,0.05),dtype="d"),
                  numpy.ndarray m_plate_num_div=numpy.array((2,2,2),dtype="i"),
                  numpy.ndarray m_plate_prop=numpy.array((500,1e7,0.3),dtype="d"),
                  numpy.ndarray m_gravity=numpy.array((0,0,-9.8),dtype="d"),
                  numpy.ndarray m_free_x=numpy.array((1,1,1),dtype="i")
                  ):
        self.my_system = _system
        self.thisptr = newChFlexPlate(self.my_system.thisptr.systemSMC,
                                     m_is3D,
                                      m_timeStep,
                                     <double*> m_plate_center.data,
                                     <double*> m_plate_dims.data,
                                     <int*> m_plate_num_div.data,
                                     <double*> m_plate_prop.data,
                                     <double*> m_gravity.data,
                                     <int*> m_free_x.data
                                     )

        self.n_coor=np.zeros((self.thisptr.num_total_nodes,3))
        self.v_coor=np.zeros((self.thisptr.num_total_nodes,3))
        self.calcNodalInfo()

    def getNumNodes(self):
        logEvent("Chrono object created with " +`self.thisptr.num_total_nodes` + " nodes")
        return self.thisptr.num_total_nodes

    def attachModel(self,model,ar):
        self.model=model
        return self

    def get_u(self):
        return 0
    def get_v(self):
        return 0
    def get_w(self):
        return 0

    def calculate_init(self):
        self.calculate()

    def SyncData(self,numpy.ndarray new_coors, numpy.ndarray new_coor_vel):
        self.thisptr.prepareData(<double*> new_coors.data,<double*> new_coor_vel.data)

    def d_N_IBM(self, numpy.ndarray x):
        import numpy as np
        cdef double dN[4]
        self.thisptr.calc_d_N_IBM(<double*> x.data,dN)
        #print(dN[0],dN[1],dN[2],dN[3])
        dist=dN[0]
        normal=np.array([dN[1],dN[2],dN[3]])
        return dist, normal

    def vel_IBM(self, numpy.ndarray x):
        import numpy as np
        cdef double vel[3]
        self.thisptr.calc_vel_IBM(<double*> x.data,vel)
        vel_np=np.array([vel[0],vel[1],vel[2]])
        return vel_np

    def calcNodalInfo(self):
        import  numpy as np

        pos=self.thisptr.calcNodalPos()
        normal=self.thisptr.calcNodalNormal()
        np_pos=np.zeros((self.thisptr.num_surface_nodes,3))
        np_normal=np.zeros((self.thisptr.num_surface_nodes,3))

        for i in range(0,self.thisptr.num_surface_nodes,1):
            np_pos[i,0]=pos[3*i+0]
            np_pos[i,1]=pos[3*i+1]
            np_pos[i,2]=pos[3*i+2]
            #print (np_pos[i,0],np_pos[i,1],np_pos[i,2])
            np_normal[i,0]=normal[3*i+0]
            np_normal[i,1]=normal[3*i+1]
            np_normal[i,2]=normal[3*i+2]
            #print (np_normal[i,0],np_normal[i,1],np_normal[i,2])
        return np_pos, np_normal


    def step(self,
             numpy.ndarray force,
             int num_nodes,
             double dt):
        self.thisptr.step(<double*> force.data,
                          num_nodes,
                          dt)

    def calculate(self,numpy.ndarray Forces, double dt):
        import  numpy as np
        from numpy.linalg import inv
        import copy
        logEvent("Calling chrono with dt " +`dt`)
        self.step(Forces,0,dt)
        # logEvent("After Step ")

        for i in range(0,self.thisptr.num_total_nodes,1):
            self.n_coor[i,0]=self.thisptr.nodal_corrd_vec[3*i+0]
            self.n_coor[i,1]=self.thisptr.nodal_corrd_vec[3*i+1]
            self.n_coor[i,2]=self.thisptr.nodal_corrd_vec[3*i+2]
            self.v_coor[i,0]=self.thisptr.nodal_corrd_vel_vec[3*i+0]
            self.v_coor[i,1]=self.thisptr.nodal_corrd_vel_vec[3*i+1]
            self.v_coor[i,2]=self.thisptr.nodal_corrd_vel_vec[3*i+2]
        return self.n_coor, self.v_coor


    def writeFrame(self):
        self.thisptr.writeThisFrame()
