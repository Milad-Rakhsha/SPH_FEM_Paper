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
// This file contains the chrono mode to set a FSI simulation.
// =============================================================================

#include "chrono/physics/ChSystemSMC.h"
#include "chrono/solver/ChSolverMINRES.h"
#include "chrono/solver/ChSolverSMC.h"
#include "chrono/fea/ChElementBrick.h"
#include "chrono/fea/ChNodeFEAxyz.h"
#include "chrono/fea/ChMesh.h"
#include "chrono_thirdparty/filesystem/path.h"
#include "outputUtil.h"

const std::string data_folder = "outputs/";
const std::string mesh_file = data_folder + "MESH";

void findNearestElements_and_NatrualCoordinates(double *pos,
                                                int &My_Chrono_Elem,
                                                ChVector<double> &Natural_coordinates,
                                                std::vector<std::vector<int>> Element_nodes,
                                                std::vector<ChVector<double>> nodal_corrd);

void find_neighbor_nodes(std::vector<std::vector<int>> &my_neighbors,
                         std::vector<std::vector<int>> elem_nodes,
                         std::vector<ChVector<double>> nodal_corrd,
                         int total_num_nodes,
                         double max_dist);

void writeMesh(std::shared_ptr<ChMesh> my_mesh, string SaveAs, std::vector<std::vector<int>> &NodeNeighborElement);
class cppChFlexPlate
{
  public:
    cppChFlexPlate(std::shared_ptr<ChSystemSMC> _system,
                   bool m_is3D,
                   double m_timeStep,
                   double *m_plate_center,
                   double *m_plate_dims,
                   int *m_plate_num_div,
                   double *m_plate_prop,
                   double *m_gravity,
                   int *m_free_x)
        : is3D(m_is3D),
          time_step(m_timeStep),
          out_step(0.02),
          free_x(m_free_x)
    {

        my_system = _system;
        const std::string rmCmd = (std::string("rm ") + data_folder + std::string("/*"));
        system(rmCmd.c_str());

        if (filesystem::create_directory(filesystem::path(data_folder.c_str())) < 0)
        {
            double a = 0;
            std::cout << "Error creating directory " << data_folder << std::endl;
            std::cin >> a;
        }
        memcpy(plate_center, m_plate_center, 3 * sizeof(double *));
        memcpy(plate_dims, m_plate_dims, 3 * sizeof(double *));
        memcpy(plate_num_div, m_plate_num_div, 3 * sizeof(int *));
        memcpy(plate_prop, m_plate_prop, 3 * sizeof(double *));
        memcpy(gravity, m_gravity, 3 * sizeof(double *));

        printf("is 3D=%d\n", m_is3D);
        // The soft-real-time cycle
        double time = 0.0;
        double out_time = 0.0;

        // The physical system: it contains all physical objects.
        // Create a mesh, that is a container for groups
        // of elements and their referenced nodes.
        my_mesh = std::make_shared<ChMesh>();
        // Geometry of the plate
        double plate_lenght_x = m_plate_dims[0];
        double plate_lenght_y = m_plate_dims[1];
        double plate_lenght_z = m_plate_dims[2]; // small thickness
        printf("plate_lenght_x=%f,plate_lenght_y=%f,plate_lenght_z=%f\n", plate_lenght_x, plate_lenght_y,
               plate_lenght_z);

        // Specification of the mesh
        int numDiv_x = m_plate_num_div[0];
        int numDiv_y = m_plate_num_div[1];
        int numDiv_z = m_plate_num_div[2];
        printf("numDiv_x=%d,numDiv_y=%d,numDiv_z=%d\n\n", numDiv_x, numDiv_y, numDiv_z);

        int N_x = numDiv_x + 1;
        int N_y = numDiv_y + 1;
        int N_z = numDiv_z + 1;
        // Number of elements in the z direction is considered as 1
        int TotalNumElements = numDiv_x * numDiv_y * numDiv_z;
        int TotalNumNodes = (numDiv_x + 1) * (numDiv_y + 1) * (numDiv_z + 1);
        num_total_nodes = TotalNumNodes;

        // For uniform mesh
        double dx = plate_lenght_x / numDiv_x;
        double dy = plate_lenght_y / numDiv_y;
        double dz = plate_lenght_z / numDiv_z;

        double density = m_plate_prop[0];
        double E = m_plate_prop[1];
        double nu = m_plate_prop[2];

        Element_nodes.resize(TotalNumElements);
        nodal_corrd_last.resize(TotalNumNodes);
        nodal_corrd.resize(TotalNumNodes);

        auto mmaterial = std::make_shared<ChContinuumElastic>();
        mmaterial->Set_RayleighDampingK(0.0);
        mmaterial->Set_RayleighDampingM(0.0);
        mmaterial->Set_density(density);
        mmaterial->Set_E(E);
        mmaterial->Set_G(E / 2 / (1 + nu));
        mmaterial->Set_v(nu);
        int numnode = 0;
        for (int k = 0; k < N_z; k++)
        {
            for (int j = 0; j < N_y; j++)
            {
                for (int i = 0; i < N_x; i++)
                {
                    double loc_x = i * dx - plate_lenght_x / 2 + m_plate_center[0];
                    double loc_y = j * dy - plate_lenght_y / 2 + m_plate_center[1];
                    double loc_z = k * dz - plate_lenght_z / 2 + m_plate_center[2];

                    ChVector<> myLoc(loc_x, loc_y, loc_z);

                    auto node = std::make_shared<ChNodeFEAxyz>(myLoc);
                    node->SetMass(0.0);

                    if (j == 0)
                        node->SetFixed(true);
                    if (j == N_y - 1)
                        node_tip = node;

                    nodal_corrd[numnode] = myLoc;
                    nodal_corrd_last[numnode] = myLoc;

                    my_mesh->AddNode(node);
                    numnode++;
                }
            }
        }

        printf("TotalNumNodes = %d\n", TotalNumNodes);
        tipnode[0] = +plate_lenght_x / 2 + m_plate_center[0];
        tipnode[1] = -plate_lenght_y / 2 + m_plate_center[1];
        tipnode[2] = -plate_lenght_z / 2 + m_plate_center[2];
        nodal_corrd_last = nodal_corrd;

        ChMatrixNM<double, 3, 1> Dims; // read element length, used in ChElementBrick
        Dims.Reset();

        Dims(0, 0) = dx;
        Dims(1, 0) = dy;
        Dims(2, 0) = dz;

        int num_elem = 0;
        for (int k = 0; k < numDiv_z; k++)
        {
            for (int j = 0; j < numDiv_y; j++)
            {
                for (int i = 0; i < numDiv_x; i++)
                {
                    auto element = std::make_shared<ChElementBrick>();
                    element->SetInertFlexVec(Dims);

                    int node0 = (i + 0) + N_x * (j + 0) + N_x * N_y * k;
                    int node1 = (i + 1) + N_x * (j + 0) + N_x * N_y * k;
                    int node2 = (i + 1) + N_x * (j + 1) + N_x * N_y * k;
                    int node3 = (i + 0) + N_x * (j + 1) + N_x * N_y * k;

                    int node4 = (i + 0) + N_x * (j + 0) + N_x * N_y * (k + 1);
                    int node5 = (i + 1) + N_x * (j + 0) + N_x * N_y * (k + 1);
                    int node6 = (i + 1) + N_x * (j + 1) + N_x * N_y * (k + 1);
                    int node7 = (i + 0) + N_x * (j + 1) + N_x * N_y * (k + 1);

                    //                    printf("\nElement %d nodes are: %d,%d,%d,%d,%d,%d,%d,%d,  ", num_elem, node0,
                    //                    node1, node2, node3,
                    //                           node4, node5, node6, node7);

                    Element_nodes[num_elem].resize(8);
                    Element_nodes[num_elem] = {node0, node1, node2, node3, node4, node5, node6, node7};

                    ChVector<double> my_center =
                        (1 / 8.0) * (nodal_corrd[Element_nodes[num_elem][0]] + nodal_corrd[Element_nodes[num_elem][1]] +
                                     nodal_corrd[Element_nodes[num_elem][2]] + nodal_corrd[Element_nodes[num_elem][3]] +
                                     nodal_corrd[Element_nodes[num_elem][4]] + nodal_corrd[Element_nodes[num_elem][5]] +
                                     nodal_corrd[Element_nodes[num_elem][6]] + nodal_corrd[Element_nodes[num_elem][7]]);

                    //                    printf("construct with center : %f,%f,%f\n", my_center.x(), my_center.y(),
                    //                    my_center.z());

                    element->SetNodes(std::dynamic_pointer_cast<ChNodeFEAxyz>(my_mesh->GetNode(node0)),
                                      std::dynamic_pointer_cast<ChNodeFEAxyz>(my_mesh->GetNode(node1)),
                                      std::dynamic_pointer_cast<ChNodeFEAxyz>(my_mesh->GetNode(node2)),
                                      std::dynamic_pointer_cast<ChNodeFEAxyz>(my_mesh->GetNode(node3)),
                                      std::dynamic_pointer_cast<ChNodeFEAxyz>(my_mesh->GetNode(node4)),
                                      std::dynamic_pointer_cast<ChNodeFEAxyz>(my_mesh->GetNode(node5)),
                                      std::dynamic_pointer_cast<ChNodeFEAxyz>(my_mesh->GetNode(node6)),
                                      std::dynamic_pointer_cast<ChNodeFEAxyz>(my_mesh->GetNode(node7)));

                    element->SetMaterial(mmaterial);
                    //    element->SetElemNum(i);       // for EAS
                    element->SetGravityOn(true);              // turn gravity on/off from within the element
                    element->SetMooneyRivlin(false);          // turn on/off Mooney Rivlin (Linear Isotropic by default)
                    ChMatrixNM<double, 9, 1> stock_alpha_EAS; //
                    stock_alpha_EAS.Reset();
                    element->SetStockAlpha(stock_alpha_EAS(0, 0), stock_alpha_EAS(1, 0), stock_alpha_EAS(2, 0),
                                           stock_alpha_EAS(3, 0), stock_alpha_EAS(4, 0), stock_alpha_EAS(5, 0),
                                           stock_alpha_EAS(6, 0), stock_alpha_EAS(7, 0), stock_alpha_EAS(8, 0));
                    my_mesh->AddElement(element);
                    num_elem++;
                }
            }
        }

        find_neighbor_nodes(Neighbors, Element_nodes, nodal_corrd, TotalNumNodes, Dims.Max());
        for (int i = 0; i < Neighbors.size(); i++)
        {
            if (Neighbors[i].size() < 6)
            {
                if (m_is3D || (!m_is3D && std::abs(nodal_corrd[i].z()) < 1e-6))
                {
                    surface_nodes.push_back(i);
                    num_surface_nodes++;
                }
            }
        }
        printf("Chrono has %d surface nodes for FS coupling", surface_nodes.size());
        Nodal_Pos = (double *)malloc(num_surface_nodes * 3 * sizeof(double *));
        Nodal_Normal = (double *)malloc(num_surface_nodes * 3 * sizeof(double *));
        nodal_corrd_vec = (double *)malloc(num_total_nodes * 3 * sizeof(double *));
        nodal_corrd_vel_vec = (double *)malloc(num_total_nodes * 3 * sizeof(double *));
        // memset(nodal_corrd_vel_vec, 0.0, sizeof(double) * num_total_nodes);

        for (int i = 0; i < num_total_nodes; i++)
        {
            nodal_corrd_vec[i * 3 + 0] = nodal_corrd[i].x();
            nodal_corrd_vec[i * 3 + 1] = nodal_corrd[i].y();
            nodal_corrd_vec[i * 3 + 2] = nodal_corrd[i].z();

            nodal_corrd_vel_vec[i * 3 + 0] = 0.0;
            nodal_corrd_vel_vec[i * 3 + 1] = 0.0;
            nodal_corrd_vel_vec[i * 3 + 2] = 0.0;
        }
        // Deactivate automatic gravity in mesh
        my_mesh->SetAutomaticGravity(false);
        my_system->Set_G_acc(ChVector<>(m_gravity[0], m_gravity[1], m_gravity[2]));

        // Remember to add the mesh to the system!
        my_system->Add(my_mesh);
        // Mark completion of system construction
        my_system->SetupInitial();

        // Perform a dynamic time integration:
        my_system->SetSolverType(ChSolver::Type::MINRES);
        auto msolver = std::static_pointer_cast<ChSolverMINRES>(my_system->GetSolver());
        msolver->SetDiagonalPreconditioning(true);
        my_system->SetMaxItersSolverSpeed(10000);
        my_system->SetTolForce(1e-08);

        my_system->SetTimestepperType(ChTimestepper::Type::EULER_IMPLICIT);

        //        my_system->SetTimestepperType(ChTimestepper::Type::HHT);
        //        auto mystepper = std::dynamic_pointer_cast<ChTimestepperHHT>(my_system->GetTimestepper());
        //        mystepper->SetAlpha(-0.2);
        //        mystepper->SetMaxiters(100);
        //        mystepper->SetAbsTolerances(1e-5);
        //        mystepper->SetMode(ChTimestepperHHT::POSITION);
        //        mystepper->SetScaling(true);
        //        mystepper->SetVerbose(true);
        if (filesystem::create_directory(filesystem::path(data_folder.c_str())) < 0)
        {
            int a;
            cout << "Error creating directory, Pres any key to continue " << data_folder << endl;
            cin >> a;
        }
        writeMesh(my_mesh, mesh_file, NodeNeighborElement);
        writeThisFrame();
    }

    bool is3D;
    double *plate_center = new double[3]; //(x,y,z)
    double *plate_dims = new double[3];   //(height,width,thickness)
    int *plate_num_div = new int[3];      //(Nx,Ny,Nz) number of division in each direction
    double *plate_prop = new double[3];   //(E,nu,G)
    double *gravity = new double[3];      //
    int *free_x;                          // For debugging purposes

    double time_step;
    double out_step;
    double time = 0.0;
    double out_time = 0.0;
    int barId;
    //    int binId;
    double tipnode[3];

    std::shared_ptr<ChSystemSMC> my_system;
    std::shared_ptr<ChMesh> my_mesh;
    std::vector<std::vector<int>> Element_nodes;       // Index of the nodes in an element
    std::vector<std::vector<int>> NodeNeighborElement; // Neighbor element of nodes in a vector
    std::vector<ChVector<double>> nodal_corrd_last;    // Position of the nodes before the update
    std::vector<ChVector<double>> nodal_corrd;         // Position of the nodes
    std::shared_ptr<ChNodeFEAxyz> node_tip;
    double *nodal_corrd_vec;
    double *nodal_corrd_vel_vec;

    std::vector<std::vector<int>> Neighbors;

    std::vector<int> surface_nodes;
    std::vector<double> surface_nodes_area;

    int num_surface_nodes = 0;
    int num_total_nodes = 0;

    int outframe = 0;
    double *Nodal_Pos;
    double *Nodal_Normal;

    // Note that forces are nodal forces in chrono. Hence size of force
    // pos are the position of the vertices in proteus
    double *step(double *forces, int num_nodes, double dt)
    {
        ChVector<double> totalForce(0);

        // // TODO: This averaged area is just a rough approximation for this problem
        // // The pressure should be integrated based on some accurate average area of each node.
        double dA = plate_dims[1] * plate_dims[2] * 2 / surface_nodes.size();

        for (int i = 0; i < surface_nodes.size(); i++)
        {
            auto node = std::dynamic_pointer_cast<ChNodeFEAxyz>(my_mesh->GetNode(surface_nodes[i]));

            ChVector<double> myForce = ChVector<double>(forces[3 * i + 0], forces[3 * i + 1], forces[3 * i + 2]);
            int numNeigh = Neighbors[i].size();
            if (numNeigh == 5)
                numNeigh--;

            node->SetForce(myForce * dA);
            totalForce += myForce * dA;
        }

        printf("totalForce to chrono dA=%f (%f, %f, %f) \n", dA, totalForce.x(), totalForce.y(), totalForce.z());

        SaveNodalCor(nodal_corrd_last);

        // my_system->DoStepDynamics(dt);

        int sync = std::round(dt / time_step);
        if (sync >= 1)
        {
            printf("%d * DoStepChronoSystem with dt= %f\n", sync, dt / sync);
            for (int t = 0; t < sync; t++)
            {
                my_system->DoStepDynamics(dt / sync);
                time += dt / sync;
            }
        }
        else
        {
            printf("DoStepChronoSystem with dt= %f\n ", dt);
            my_system->DoStepDynamics(dt);
        }

        SaveNodalCor(nodal_corrd);

        for (int i = 0; i < nodal_corrd.size(); i++)
        {
            ChVector<> pos = std::dynamic_pointer_cast<ChNodeFEAxyz>(my_mesh->GetNode(i))->GetPos();
            nodal_corrd_vec[i * 3 + 0] = pos.x();
            nodal_corrd_vec[i * 3 + 1] = pos.y();
            nodal_corrd_vec[i * 3 + 2] = pos.z();
            ChVector<> vel = std::dynamic_pointer_cast<ChNodeFEAxyz>(my_mesh->GetNode(i))->GetPos_dt();
            nodal_corrd_vel_vec[i * 3 + 0] = vel.x();
            nodal_corrd_vel_vec[i * 3 + 1] = vel.y();
            nodal_corrd_vel_vec[i * 3 + 2] = vel.z();
        }

        printf("chrono time : %f\n", my_system->GetChTime());
        return nodal_corrd_vec;
    }

    void writeThisFrame()
    {
        char SaveAsBuffer[256]; // The filename buffer.
        snprintf(SaveAsBuffer, sizeof(char) * 256, (data_folder + "flex_body.%d.vtk").c_str(), outframe);
        writeFrame(my_mesh, SaveAsBuffer, mesh_file, NodeNeighborElement);
        outframe++;
    }
    // This function must be called to update the std::vector representation of the nodal coordinates
    // This is important when a rank broadcast information to others, and processor start calculating hx,hy,hz
    // This must be called before proceeding to hx,hy,hz
    void prepareData(double *nodal_corrd_v, double *nodal_corrd_vel_v)
    {
        for (int i = 0; i < nodal_corrd.size(); i++)
        {
            nodal_corrd[i].x() = nodal_corrd_v[i * 3 + 0];
            nodal_corrd[i].y() = nodal_corrd_v[i * 3 + 1];
            nodal_corrd[i].z() = nodal_corrd_v[i * 3 + 2];
            nodal_corrd_vel_vec[i * 3 + 0] = nodal_corrd_vel_v[i * 3 + 0];
            nodal_corrd_vel_vec[i * 3 + 1] = nodal_corrd_vel_v[i * 3 + 1];
            nodal_corrd_vel_vec[i * 3 + 2] = nodal_corrd_vel_v[i * 3 + 2];
        }
    }

    double *calcNodalPos()
    {
        for (int i = 0; i < surface_nodes.size(); i++)
        {
            ChVector<> pos = nodal_corrd[surface_nodes[i]];
            Nodal_Pos[i * 3 + 0] = pos.x();
            Nodal_Pos[i * 3 + 1] = pos.y();
            Nodal_Pos[i * 3 + 2] = pos.z();
        }
        return Nodal_Pos;
    }

    double *calcNodalNormal()
    {
        for (int i = 0; i < surface_nodes.size(); i++)
        {
            int mynode = surface_nodes[i];
            ChVector<> normal(0);
            for (int j = 0; j < Neighbors[mynode].size(); j++)
            {
                int nei = Neighbors[mynode][j];
                double L = (nodal_corrd[mynode] - nodal_corrd[nei]).Length();
                normal += (nodal_corrd[mynode] - nodal_corrd[nei]) / L;
            }
            normal.Normalize();
            Nodal_Normal[i * 3 + 0] = normal.x();
            Nodal_Normal[i * 3 + 1] = normal.y();
            Nodal_Normal[i * 3 + 2] = normal.z();
        }
        return Nodal_Normal;
    }

    void calcPos(ChMatrix<> N, int Chrono_elem, ChVector<> &new_pos, std::vector<ChVector<double>> nod_cor)
    {
        ChVector<double> p0 = nod_cor[Element_nodes[Chrono_elem][0]];
        ChVector<double> p1 = nod_cor[Element_nodes[Chrono_elem][1]];
        ChVector<double> p2 = nod_cor[Element_nodes[Chrono_elem][2]];
        ChVector<double> p3 = nod_cor[Element_nodes[Chrono_elem][3]];
        ChVector<double> p4 = nod_cor[Element_nodes[Chrono_elem][4]];
        ChVector<double> p5 = nod_cor[Element_nodes[Chrono_elem][5]];
        ChVector<double> p6 = nod_cor[Element_nodes[Chrono_elem][6]];
        ChVector<double> p7 = nod_cor[Element_nodes[Chrono_elem][7]];
        new_pos = N(0) * p0 + N(1) * p1 + N(2) * p2 + N(3) * p3 + N(4) * p4 + N(5) * p5 + N(6) * p6 + N(7) * p7;
    }

    void calcVel(ChMatrix<> N, int Chrono_elem, ChVector<> &new_vel)
    {
        std::vector<int> thisElem = Element_nodes[Chrono_elem];
        int node = thisElem[0];
        ChVector<> v0(nodal_corrd_vel_vec[3 * node + 0], nodal_corrd_vel_vec[3 * node + 1], nodal_corrd_vel_vec[3 * node + 2]);
        node = thisElem[1];
        ChVector<> v1(nodal_corrd_vel_vec[3 * node + 0], nodal_corrd_vel_vec[3 * node + 1], nodal_corrd_vel_vec[3 * node + 2]);
        node = thisElem[2];
        ChVector<> v2(nodal_corrd_vel_vec[3 * node + 0], nodal_corrd_vel_vec[3 * node + 1], nodal_corrd_vel_vec[3 * node + 2]);
        node = thisElem[3];
        ChVector<> v3(nodal_corrd_vel_vec[3 * node + 0], nodal_corrd_vel_vec[3 * node + 1], nodal_corrd_vel_vec[3 * node + 2]);
        node = thisElem[4];
        ChVector<> v4(nodal_corrd_vel_vec[3 * node + 0], nodal_corrd_vel_vec[3 * node + 1], nodal_corrd_vel_vec[3 * node + 2]);
        node = thisElem[5];
        ChVector<> v5(nodal_corrd_vel_vec[3 * node + 0], nodal_corrd_vel_vec[3 * node + 1], nodal_corrd_vel_vec[3 * node + 2]);
        node = thisElem[6];
        ChVector<> v6(nodal_corrd_vel_vec[3 * node + 0], nodal_corrd_vel_vec[3 * node + 1], nodal_corrd_vel_vec[3 * node + 2]);
        node = thisElem[7];
        ChVector<> v7(nodal_corrd_vel_vec[3 * node + 0], nodal_corrd_vel_vec[3 * node + 1], nodal_corrd_vel_vec[3 * node + 2]);
        new_vel = N(0) * v0 + N(1) * v1 + N(2) * v2 + N(3) * v3 + N(4) * v4 + N(5) * v5 + N(6) * v6 + N(7) * v7;
    }
    void SaveNodalCor(std::vector<ChVector<double>> &nod_cor)
    {
        for (int i = 0; i < my_mesh->GetNnodes(); i++)
        {
            nod_cor[i] = std::dynamic_pointer_cast<ChNodeFEAxyz>(my_mesh->GetNode(i))->GetPos();
        }
    }

    void calc_d_N_IBM(double *x, double *output)
    {
        double d = 1e6;
        ChVector<> p(x[0], x[1], x[2]);
        ChVector<> N(0.0);
        ChVector<> d_vector(0.0);
        for (int i = 0; i < surface_nodes.size(); i++)
        {
            int nodenum = surface_nodes[i];
            double L = (p - nodal_corrd[nodenum]).Length();
            if (L < d)
            {
                d = L;
                N = ChVector<>(Nodal_Normal[3 * i + 0], Nodal_Normal[3 * i + 1], Nodal_Normal[3 * i + 2]);
                d_vector = (p - nodal_corrd[nodenum]);
                // output[0] = d_vector.Dot(N) / N.Length();
                if (d_vector.Dot(N) / N.Length() > 0)
                    output[0] = d_vector.Length();
                else
                    output[0] = -d_vector.Length();
                d_vector /= d_vector.Length();
                output[1] = d_vector.x(); //Nodal_Normal[0];
                output[2] = d_vector.y(); //Nodal_Normal[1];
                output[3] = d_vector.z(); //Nodal_Normal[2];
            }
        }
    }

    void calc_vel_IBM(double *x, double *output)
    {
        double d = 1e6;
        double Characteristic_L = std::abs(nodal_corrd[0].x() - nodal_corrd[1].x());
        ChVector<>
            p(x[0], x[1], x[2]);
        ChVector<> d_vector(0.0);

        int closestnode = -1;
        for (int i = 0; i < num_total_nodes; i++)
        {
            double L = (p - nodal_corrd[i]).Length();
            if (L < d && L < 2 * Characteristic_L)
            {
                d = L;
                closestnode = i;
            }
        }
        if (closestnode != -1)
        {
            int closestElem = 0;
            ChVector<> Natural_coordinates(0);
            findNearestElements_and_NatrualCoordinates(x, closestElem, Natural_coordinates, Element_nodes, nodal_corrd);
            auto elem = std::dynamic_pointer_cast<ChElementBrick>(my_mesh->GetElement(closestElem));
            ChMatrixNM<double, 1, 8> N;
            elem->ShapeFunctions(N, Natural_coordinates.x(), Natural_coordinates.y(), Natural_coordinates.z());
            // printf("my elem=%d\t",closestElem);// Natural_coordinates=%f,%f,%f",closestElem,Natural_coordinates.x(),Natural_coordinates.y(),Natural_coordinates.z());
            ChVector<> out_vel(0);
            calcVel(N, closestElem, out_vel);
            output[0] = out_vel.x();
            output[1] = out_vel.y();
            output[2] = out_vel.z();
        }
    }
};

/////////////////////////////////////////////////////////
/////// End of the class cppChFlexPlate//////////////////
/////////////////////////////////////////////////////////

cppChFlexPlate *newChFlexPlate(std::shared_ptr<ChSystemSMC> _system,
                               bool m_is3D,
                               double m_timeStep,
                               double *m_plate_center,
                               double *m_plate_dims,
                               int *m_plate_num_div,
                               double *m_plate_prop,
                               double *m_gravity,
                               int *m_free_x)
{
    return new cppChFlexPlate(_system, m_is3D, m_timeStep, m_plate_center, m_plate_dims, m_plate_num_div, m_plate_prop,
                              m_gravity, m_free_x);
}

// This function gives a vector of points, searches in the chrono elements and return a vector of "the nearest
// element of each point"
void findNearestElements_and_NatrualCoordinates(double *pos,
                                                int &My_Chrono_Elem,
                                                ChVector<double> &Natural_coordinates,
                                                std::vector<std::vector<int>> elem_nodes,
                                                std::vector<ChVector<double>> nod_cor)
{
    std::vector<ChVector<double>> Elem_centers;
    for (int i = 0; i < elem_nodes.size(); i++)
    {
        ChVector<double> my_center =
            (1 / 8.0) * (nod_cor[elem_nodes[i][0]] + nod_cor[elem_nodes[i][1]] + nod_cor[elem_nodes[i][2]] +
                         nod_cor[elem_nodes[i][3]] + nod_cor[elem_nodes[i][4]] + nod_cor[elem_nodes[i][5]] +
                         nod_cor[elem_nodes[i][6]] + nod_cor[elem_nodes[i][7]]);

        Elem_centers.push_back(my_center);
    }

    double min = 1e20;
    int nearest_elem = 0;
    ChVector<double> mypos(pos[0], pos[1], pos[2]);
    for (int i = 0; i < Elem_centers.size(); i++)
    {
        double dist = (mypos - Elem_centers[i]).Length();
        if (dist <= min)
        {
            min = dist;
            nearest_elem = i;
        }
    }
    My_Chrono_Elem = nearest_elem;

    // d3 is in the global coordinate frame
    ChVector<double> d3 = mypos - Elem_centers[nearest_elem];
    // need to find d3 in the element local coordinate
    ChVector<double> x_dir = (+nod_cor[elem_nodes[nearest_elem][1]] + nod_cor[elem_nodes[nearest_elem][2]] +
                              +nod_cor[elem_nodes[nearest_elem][5]] + nod_cor[elem_nodes[nearest_elem][6]] +
                              -nod_cor[elem_nodes[nearest_elem][0]] - nod_cor[elem_nodes[nearest_elem][3]] +
                              -nod_cor[elem_nodes[nearest_elem][4]] - nod_cor[elem_nodes[nearest_elem][7]]) /
                             8.0;

    ChVector<double> y_dir = (+nod_cor[elem_nodes[nearest_elem][2]] + nod_cor[elem_nodes[nearest_elem][3]] +
                              +nod_cor[elem_nodes[nearest_elem][6]] + nod_cor[elem_nodes[nearest_elem][7]] +
                              -nod_cor[elem_nodes[nearest_elem][0]] - nod_cor[elem_nodes[nearest_elem][1]] +
                              -nod_cor[elem_nodes[nearest_elem][4]] - nod_cor[elem_nodes[nearest_elem][5]]) /
                             8.0;

    ChVector<double> z_dir = (+nod_cor[elem_nodes[nearest_elem][4]] + nod_cor[elem_nodes[nearest_elem][5]] +
                              +nod_cor[elem_nodes[nearest_elem][6]] + nod_cor[elem_nodes[nearest_elem][7]] +
                              -nod_cor[elem_nodes[nearest_elem][0]] - nod_cor[elem_nodes[nearest_elem][1]] +
                              -nod_cor[elem_nodes[nearest_elem][2]] - nod_cor[elem_nodes[nearest_elem][3]]) /
                             8.0;

    ChVector<double> myNatural_coordinate;
    double x_l = x_dir.Length();
    x_dir.Normalize();
    double y_l = y_dir.Length();
    y_dir.Normalize();
    double z_l = z_dir.Length();
    z_dir.Normalize();

    //    printf("\n x_d=%f,%f,%f, y_d=%f,%f,%f, z_d=%f,%f,%f\n", x_dir.x(), x_dir.y(), x_dir.z(), y_dir.x(),
    //    y_dir.y(),
    //           y_dir.z(), z_dir.x(), z_dir.y(), z_dir.z());
    myNatural_coordinate.x() = d3.Dot(x_dir) / x_l;
    myNatural_coordinate.y() = d3.Dot(y_dir) / y_l;
    myNatural_coordinate.z() = d3.Dot(z_dir) / z_l;

    //    printf("d3=%f,%f,%f, n Cor= %f,%f,%f\n", d3.x(), d3.y(), d3.z(), myNatural_coordinate.x(),
    //    myNatural_coordinate.y(),
    //           myNatural_coordinate.z());

    Natural_coordinates = myNatural_coordinate;
};

bool isNeightbor(int i, int j)
{
    if (i == 0 && (j == 1 || j == 3 || j == 4))
        return true;
    if (i == 1 && (j == 0 || j == 2 || j == 5))
        return true;
    if (i == 2 && (j == 1 || j == 3 || j == 6))
        return true;
    if (i == 3 && (j == 0 || j == 2 || j == 7))
        return true;
    if (i == 4 && (j == 5 || j == 7 || j == 0))
        return true;
    if (i == 5 && (j == 6 || j == 4 || j == 1))
        return true;
    if (i == 6 && (j == 7 || j == 5 || j == 2))
        return true;
    if (i == 7 && (j == 4 || j == 6 || j == 3))
        return true;

    return false;
};
void find_neighbor_nodes(std::vector<std::vector<int>> &my_neighbors,
                         std::vector<std::vector<int>> elem_nodes,
                         std::vector<ChVector<double>> nodal_corrd,
                         int total_num_nodes,
                         double max_dist)
{
    // Iterate through the elements nodes and figure out the neighbor nodes of each node
    my_neighbors.resize(total_num_nodes);

    for (int i = 0; i < my_neighbors.size(); i++)
    {
        printf("node %d ", i);
        /////////////////////

        for (int j = 0; j < elem_nodes.size(); j++)
        {
            // In this element's nodes
            if (std::find(elem_nodes[j].begin(), elem_nodes[j].end(), i) != elem_nodes[j].end())
            {
                vector<int>::iterator it = find(elem_nodes[j].begin(), elem_nodes[j].end(), i);
                int pos = distance(elem_nodes[j].begin(), it);
                /////////////
                for (int k = 0; k < elem_nodes[j].size(); k++)
                {
                    if (!(std::find(my_neighbors[i].begin(), my_neighbors[i].end(), elem_nodes[j][k]) !=
                          my_neighbors[i].end()))
                    {
                        ///////////
                        // Check if the two nodes are neightbors
                        if (isNeightbor(pos, k))
                        {
                            my_neighbors[i].push_back(elem_nodes[j][k]);
                            printf("%d, ", elem_nodes[j][k]);
                        }
                        ////////////////////////////////////////
                    }
                }
                ///////////////////////
            }
            /////////////////////
        }
        printf("\n");
    }

    ////////////
};
