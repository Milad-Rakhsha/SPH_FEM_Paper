// =============================================================================
// PROJECT CHRONO - http://projectchrono.org
//
// Copyright (c) 2014 projectchrono.org
// All right reserved.
//
// Use of this source code is governed by a BSD-style license that can be found
// in the LICENSE file at the top level of the distribution and at
// http://projectchrono.org/license-chrono.txt.
//
// =============================================================================
// Authors: Milad Rakhsha
// =============================================================================
//
// This file contains the chrono model to set up a FSI simulation.
// =============================================================================
#include <cmath>
#include <cstdio>
#include <vector>

#include "chrono/assets/ChTriangleMeshShape.h"
#include "chrono/geometry/ChTriangleMeshConnected.h"
#include "chrono/physics/ChBodyEasy.h"
#include "chrono/physics/ChParticlesClones.h"
#include "chrono/physics/ChSystemNSC.h"
#include "chrono/utils/ChUtilsCreators.h"
#include "chrono/utils/ChUtilsGenerators.h"
#include "chrono/utils/ChUtilsGeometry.h"
#include "chrono/utils/ChUtilsInputOutput.h"
#include "chrono_thirdparty/filesystem/path.h"
#include "chrono_thirdparty/filesystem/resolver.h"


using namespace chrono;
using namespace chrono::collision;
std::string out_dir = "outputs/";
double time_step = 1e-3;

class cppMBDModel {
 private:
  double *container_center = new double[3];  //(x,y,z)
  double *container_dims = new double[3];    //(height,width,thickness)
  double wallThickness = 0.1;

  double *particles_center = new double[3];  //(x,y,z)
  double *particles_dims = new double[3];    //(height,width,thickness)
  double diameter = 0.24;

  double *gravity = new double[3];

  ChSystemNSC my_system;

  double time_step;
  int outframe = 0;
  double time = 0.0;

 public:
  // chrono::geometry::ChTriangleMeshConnected mmeshbox;
  std::shared_ptr<geometry::ChTriangleMeshConnected> mmeshbox = std::make_shared<geometry::ChTriangleMeshConnected>();

  std::shared_ptr<ChBody> mbody;
  int objPoints = 0;
  ChMatrix33<> ROT;
  ChVector<> POS;
  ChVector<> ROT_dt;
  ChVector<> POS_dt;
  std::vector<ChVector<double>> nc;
  std::vector<ChVector<double>> nn;

 public:
  cppMBDModel(double m_timeStep, double *m_container_center,
              double *m_container_dims, double *m_particles_center,
              double *m_particles_dims, double m_particles_density,
              double m_particles_diameter, double *m_gravity)
      : time_step(m_timeStep) {
    memcpy(container_center, m_container_center, 3 * sizeof(double *));
    memcpy(container_dims, m_container_dims, 3 * sizeof(double *));
    memcpy(particles_center, m_particles_center, 3 * sizeof(double *));
    memcpy(particles_dims, m_particles_dims, 3 * sizeof(double *));
    memcpy(gravity, m_gravity, 3 * sizeof(double *));

    printf("container_dims -first= %f,%f,%f\n", container_dims[0],
           container_dims[1], container_dims[2]);

    my_system.Set_G_acc(ChVector<>(gravity[0], gravity[1], gravity[2]));
    std::cout << "Set_G_acc " << gravity[0] << " " << gravity[1] << " "
              << gravity[2] << " " << std::endl;

    const std::string rmCmd =
        (std::string("rm ") + out_dir + std::string("/*"));
    system(rmCmd.c_str());



      if (!filesystem::create_directory(filesystem::path(out_dir))) {
    		std::cout << "Error creating directory " << out_dir << std::endl;
    	}

    ChVector<double> origin(0, 0, 0.0);
    ChVector<double> parCenter(particles_center[0], particles_center[1],
                               particles_center[2]);
    ChVector<double> parHalfDim(particles_dims[0] / 2, particles_dims[1] / 2,
                                particles_dims[2] / 2);

    mbody = std::make_shared<ChBody>();
    mbody->SetPos(parCenter);
    double Vol = 3.1415 * m_particles_dims[0] * m_particles_dims[0] / 4.0 * m_particles_dims[1];
    mbody->SetMass(m_particles_density * Vol);


    mbody->SetBodyFixed(false);
    ROT = ChMatrix33<>(1.0);
    POS = parCenter;
    ROT_dt = ChVector<>(0.0, 0.0, 0.0);
    POS_dt = ChVector<>(0.0, 0.0, 0.0);

    mmeshbox->LoadWavefrontMesh("cylinder.obj", true, false);
    nc = mmeshbox->m_vertices;
    nn = mmeshbox->m_normals;

    auto mat = std::make_shared<ChMaterialSurfaceNSC>();
    mat->SetFriction(0.1f);

    objPoints = nc.size();
    printf("an obj file imported with %d vertices\n", objPoints);
    mbody->GetCollisionModel()->ClearModel();
    mbody->GetCollisionModel()->AddTriangleMesh(
        mat, mmeshbox, false, false, parCenter, Q_from_AngAxis(0.0, VECT_X), 0.005);
    // mbody->GetMaterialSurfaceNSC()->SetFriction(0.01f);
    mbody->GetCollisionModel()->BuildModel();
    mbody->SetCollide(true);
    my_system.Add(mbody);


    // The inconsistency here is to resolve the chrono default container
    // creation. Chrono container is at z=dimZ/2 not 0
    ChVector<> center(container_center[0], container_center[1],
                      container_center[2] - container_dims[2] / 2);
    ChVector<> boxDim(container_dims[0] / 2, container_dims[1] / 2,
                      container_dims[2] / 2);
    printf("container_dims = %f,%f,%f\n", container_dims[0], container_dims[1],
           container_dims[2]);
    printf("container_center = %f,%f,%f\n", container_center[0],
           container_center[1], container_center[2]);
    printf("boxDim = %f,%f,%f\n", boxDim.x(), boxDim.y(), boxDim.z());
    printf("center = %f,%f,%f\n", center.x(), center.y(), center.z());
    printf("POS = %f,%f,%f\n", POS.x(), POS.y(), POS.z());
    printf("POS_dt = %f,%f,%f\n", POS_dt.x(), POS_dt.y(), POS_dt.z());
    printf("ROT_dt = %f,%f,%f\n", ROT_dt.x(), ROT_dt.y(), ROT_dt.z());

    utils::CreateBoxContainer(&my_system, 0, mat, boxDim, wallThickness, center,
                              Q_from_AngAxis(0, VECT_Y), true, false, true,
                              true);
    writeThisFrame();
  }

  void step(double *forces, double dt) {
    ChVector<double> myForce =
        ChVector<double>(forces[0], forces[1], forces[2]);
    mbody->Empty_forces_accumulators();
    ChVector<> pos = mbody->GetPos();
    ChVector<> W = mbody->GetMass() * my_system.Get_G_acc();
    mbody->Accumulate_force(myForce, pos, false);
    printf("Add force=%f,%f,%f, mg=%f,%f,%f to the body\n", myForce.x(),
           myForce.y(), myForce.z(), W.x(), W.y(), W.z());

    int sync = std::round(dt / time_step);

    if (sync >= 1) {
      printf("%d * DoStepChronoSystem with dt= %f\n", sync, dt / sync);
      for (int i = 0; i < sync; i++) {
        my_system.DoStepDynamics(dt / sync);
        time += dt / sync;
      }
    } else {
      printf("--DoStepChronoSystem with dt= %f\n ", dt);
      my_system.DoStepDynamics(dt);
      time += dt;
    }
    mbody->Empty_forces_accumulators();

    ROT = mbody->GetRot();
    POS = mbody->GetPos();
    ROT_dt = mbody->GetWvel_loc();
    POS_dt = mbody->GetPos_dt();

    printf("===================================================\n");
    printf("chrono time : %f\n", my_system.GetChTime());
    printf("===================================================\n");
  }

  void calc_d_N_IBM(double *pos, double *out) {
    // printf("find = %f,%f,%f", pos[0], pos[1], pos[2]);

    ChVector<> p(pos[0], pos[1], pos[2]);
    ChVector<> N(0.0);
    // double output[4];
    double L = 1.0e6;
    // printf("nc.size() =%d\n", nc.size());

    for (int i = 0; i < nc.size(); i++) {
      ChVector<> d_vector = (p - (POS + ROT * nc[i]));
      double d = d_vector.Length();
      if (d < L) {
        L = d;
        N = ROT *
            nc[i];  // This is not normal vec, just the local position vector
        if (d_vector.Dot(N) > 0)
          out[0] = d_vector.Length();
        else
          out[0] = -d_vector.Length();

        d_vector /= d_vector.Length();
        out[1] = d_vector.x();
        out[2] = d_vector.y();
        out[3] = d_vector.z();
        // printf("d_vector = %f,%f,%f,%f\n", out[0], out[1], out[2], out[3]);
      }
    }
    // return output;
  }

  void calc_v_IBM(double *pos, double *out) {
    ChVector<> p(pos[0], pos[1], pos[2]);
    ChVector<> vel = POS_dt + ROT_dt % (p - POS);
    out[0] = vel.x();
    out[1] = vel.y();
    out[2] = vel.z();

    // return out;
  }

  void writeThisFrame() {
    const std::string bodies = out_dir + std::string("bodies") +
                               std::to_string(outframe) + std::string(".csv");
    utils::WriteBodies(&my_system, bodies, false, true, ",");
    outframe++;
  }
};

cppMBDModel *newMBDModel(double m_timeStep, double *m_container_center,
                         double *m_container_dims, double *m_particles_center,
                         double *m_particles_dims, double density,
                         double diameter, double *m_gravity) {
  return new cppMBDModel(m_timeStep, m_container_center, m_container_dims,
                         m_particles_center, m_particles_dims, density,
                         diameter, m_gravity);
};

int findNearestParticle(double *pos) {}
