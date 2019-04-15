import numpy
cimport numpy
from proteus import AuxiliaryVariables, Archiver
from proteus.Profiling import  logEvent
from cython.view cimport array as cvarray
import numpy as np
# Import the C-level symbols of numpy
cimport numpy as np
from libcpp cimport bool

cdef extern from "Chrono.h":
    cdef cppclass cppMBDModel:
        void step(double* forces, double dt)
        void calc_d_N_IBM(double* x, double* out)
        void calc_v_IBM(double* x, double* out)
        void writeThisFrame()
        int objPoints


    cppMBDModel* newMBDModel(double m_timeStep,
                             double* m_container_center,
                             double* m_container_dims,
                             double* m_particles_center,
                             double* m_particles_dims,
                             double m_particles_density,
                             double m_particles_diameter,
                             double* m_gravity)

cdef class MBDModel:
    cdef cppMBDModel* thisptr
    cdef object model
    cdef object numnodes

    def __cinit__(self,
                  double m_timeStep=5e-3,
                  numpy.ndarray m_container_center=numpy.array((0.0,0.0,0.0),dtype="d"),
                  numpy.ndarray m_container_dims=numpy.array((1.0,1.0,1.0),dtype="d"),
                  numpy.ndarray m_particles_center=numpy.array((0,0,0.5),dtype="d"),
                  numpy.ndarray m_particles_dims=numpy.array((0.1,0.1,0.1),dtype="d"),
                  double m_particles_density=1000,
                  double m_particles_diameter=0.05,
                  numpy.ndarray m_gravity=numpy.array((0.0,0.0,0.0),dtype="d")
                  ):

        self.thisptr =  newMBDModel(m_timeStep,
                                    <double*> m_container_center.data,
                                    <double*> m_container_dims.data,
                                    <double*> m_particles_center.data,
                                    <double*> m_particles_dims.data,
                                    m_particles_density,
                                    m_particles_diameter,
                                    <double*> m_gravity.data
                                     )

    def attachModel(self,model,ar):
        self.model=model
        return self

    def calculate_init(self):
        self.calculate()

    def get_num_surface_points(self):
        return self.thisptr.objPoints

    def calc_d_N_IBM(self, numpy.ndarray x):
        import numpy as np
        cdef double dN[4]
        self.thisptr.calc_d_N_IBM(<double*> x.data,dN)
        dist=dN[0]
        normal=np.array([dN[1],dN[2],dN[3]])
        return dist, normal


    def calc_v_IBM(self, numpy.ndarray pos):

        import numpy as np
        cdef double v[3]
        self.thisptr.calc_v_IBM(<double*> pos.data,v)
        vn=np.array([v[0],v[1],v[2]])
        return vn

    def step(self,
             numpy.ndarray force,
             double dt):

        self.thisptr.step(<double*> force.data, dt)


    def calculate(self,numpy.ndarray Forces, double dt):
        import  numpy as np
        from numpy.linalg import inv
        import copy
        logEvent("Calling chrono with dt " +`dt`)
        self.step(Forces,dt)


    def writeFrame(self):
        self.thisptr.writeThisFrame()
