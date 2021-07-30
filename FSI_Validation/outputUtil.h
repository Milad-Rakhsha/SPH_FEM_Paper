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
// This file contains the utilities to write chrono files as .vtk for paraview.
// =============================================================================

#ifndef OUTPUTUTIL_H
#define OUTPUTUTIL_H

#include "chrono/core/ChMatrix.h"
#include "chrono/physics/ChSystemSMC.h"
#include "chrono/physics/ChSystemNSC.h"
#include "chrono/fea/ChMesh.h"
#include "chrono/fea/ChElementCableANCF.h"
#include "chrono/fea/ChElementShellANCF.h"
#include "chrono/fea/ChElementBrick.h"
#include "chrono/fea/ChNodeFEAxyz.h"

using namespace chrono;
using namespace chrono::fea;
using namespace std;

void writeMesh(std::shared_ptr<ChMesh> my_mesh, string SaveAs, std::vector<std::vector<int>> &NodeNeighborElement)
{
    std::ofstream MESH; // output file stream
    MESH.open(SaveAs, std::ios::out);
    MESH.precision(7);
    MESH << std::scientific;
    NodeNeighborElement.resize(my_mesh->GetNnodes());
    std::vector<std::vector<int>> CableElemNodes;
    std::vector<std::vector<int>> ShellElemNodes;
    std::vector<std::vector<int>> BrickElemNodes;

    std::vector<std::shared_ptr<ChNodeFEAbase>> myvector;
    myvector.resize(my_mesh->GetNnodes());

    for (int i = 0; i < my_mesh->GetNnodes(); i++)
    {
        myvector[i] = std::dynamic_pointer_cast<ChNodeFEAbase>(my_mesh->GetNode(i));
    }

    int numCables = 0;
    int numShells = 0;
    int numBricks = 0;
    for (int iele = 0; iele < my_mesh->GetNelements(); iele++)
    {
        if (auto element = std::dynamic_pointer_cast<ChElementCableANCF>(my_mesh->GetElement(iele)))
            numCables++;
        if (auto element = std::dynamic_pointer_cast<ChElementShellANCF>(my_mesh->GetElement(iele)))
            numShells++;
        if (auto element = std::dynamic_pointer_cast<ChElementBrick>(my_mesh->GetElement(iele)))
            numBricks++;
    }
    MESH << "\nCELLS " << my_mesh->GetNelements() << " " << (unsigned int)(numCables * 3 + numShells * 5 + numBricks * 9) << "\n";

    for (int iele = 0; iele < my_mesh->GetNelements(); iele++)
    {
        std::vector<int> mynodes;

        if (auto element = std::dynamic_pointer_cast<ChElementCableANCF>(my_mesh->GetElement(iele)))
        {
            mynodes.resize(2);
            MESH << "2 ";
            int nodeOrder[] = {0, 1};
            mynodes[0] = element->GetNodeN(nodeOrder[0])->GetIndex();
            mynodes[1] = element->GetNodeN(nodeOrder[1])->GetIndex();
            CableElemNodes.push_back(mynodes);
            for (int myNodeN = 0; myNodeN < mynodes.size(); myNodeN++)
            {
                auto nodeA = (element->GetNodeN(nodeOrder[myNodeN]));
                std::vector<std::shared_ptr<ChNodeFEAbase>>::iterator it;
                it = find(myvector.begin(), myvector.end(), nodeA);
                if (it == myvector.end())
                {
                    // name not in vector
                }
                else
                {
                    auto index = std::distance(myvector.begin(), it);
                    MESH << (unsigned int)index << " ";
                    NodeNeighborElement[index].push_back(iele);
                }
            }
            MESH << "\n";
        }
        else if (auto element = std::dynamic_pointer_cast<ChElementShellANCF>(my_mesh->GetElement(iele)))
        {
            mynodes.resize(4);
            MESH << "4 ";
            int nodeOrder[] = {0, 1, 2, 3};
            mynodes[0] = element->GetNodeN(nodeOrder[0])->GetIndex();
            mynodes[1] = element->GetNodeN(nodeOrder[1])->GetIndex();
            mynodes[2] = element->GetNodeN(nodeOrder[2])->GetIndex();
            mynodes[3] = element->GetNodeN(nodeOrder[3])->GetIndex();
            ShellElemNodes.push_back(mynodes);
            for (int myNodeN = 0; myNodeN < mynodes.size(); myNodeN++)
            {
                auto nodeA = (element->GetNodeN(nodeOrder[myNodeN]));
                std::vector<std::shared_ptr<ChNodeFEAbase>>::iterator it;
                it = find(myvector.begin(), myvector.end(), nodeA);
                if (it == myvector.end())
                {
                    // name not in vector
                }
                else
                {
                    auto index = std::distance(myvector.begin(), it);
                    MESH << (unsigned int)index << " ";
                    NodeNeighborElement[index].push_back(iele);
                }
            }
            MESH << "\n";
        }
        else if (auto element = std::dynamic_pointer_cast<ChElementBrick>(my_mesh->GetElement(iele)))
        {
            mynodes.resize(8);
            MESH << "8 ";
            int nodeOrder[] = {0, 1, 2, 3, 4, 5, 6, 7};
            mynodes[0] = element->GetNodeN(nodeOrder[0])->GetIndex();
            mynodes[1] = element->GetNodeN(nodeOrder[1])->GetIndex();
            mynodes[2] = element->GetNodeN(nodeOrder[2])->GetIndex();
            mynodes[3] = element->GetNodeN(nodeOrder[3])->GetIndex();
            mynodes[4] = element->GetNodeN(nodeOrder[4])->GetIndex();
            mynodes[5] = element->GetNodeN(nodeOrder[5])->GetIndex();
            mynodes[6] = element->GetNodeN(nodeOrder[6])->GetIndex();

            mynodes[7] = element->GetNodeN(nodeOrder[7])->GetIndex();
            BrickElemNodes.push_back(mynodes);
            for (int myNodeN = 0; myNodeN < mynodes.size(); myNodeN++)
            {
                auto nodeA = (element->GetNodeN(nodeOrder[myNodeN]));
                std::vector<std::shared_ptr<ChNodeFEAbase>>::iterator it;
                it = find(myvector.begin(), myvector.end(), nodeA);
                if (it == myvector.end())
                {
                    // name not in vector
                }
                else
                {
                    auto index = std::distance(myvector.begin(), it);
                    MESH << (unsigned int)index << " ";
                    NodeNeighborElement[index].push_back(iele);
                }
            }
            MESH << "\n";
        }
    }

    MESH << "\nCELL_TYPES " << my_mesh->GetNelements() << "\n";

    for (int iele = 0; iele < my_mesh->GetNelements(); iele++)
    {
        if (auto element = std::dynamic_pointer_cast<ChElementCableANCF>(my_mesh->GetElement(iele)))
            MESH << "3\n";
        else if (auto element = std::dynamic_pointer_cast<ChElementShellANCF>(my_mesh->GetElement(iele)))
            MESH << "9\n";
        else if (auto element = std::dynamic_pointer_cast<ChElementBrick>(my_mesh->GetElement(iele)))
            MESH << "12\n";
    }

    // MESH.close();
    //MESH.write_to_file(SaveAs);
};

void writeFrame(std::shared_ptr<ChMesh> my_mesh,
                char SaveAsBuffer[256],
                std::string MeshFileBuffer,
                std::vector<std::vector<int>> NodeNeighborElement)
{

    std::ofstream output;
    std::string SaveAsBuffer_string(SaveAsBuffer);
    SaveAsBuffer_string.erase(SaveAsBuffer_string.length() - 4, 4);
    cout << SaveAsBuffer_string << endl;
    snprintf(SaveAsBuffer, sizeof(char) * 256, ("%s"), (SaveAsBuffer_string + ".vtk").c_str());
    output.open(SaveAsBuffer, std::ios::app);

    output << "# vtk DataFile Version 2.0" << std::endl;
    output << "Unstructured Grid Example" << std::endl;
    output << "ASCII" << std::endl;
    output << "DATASET UNSTRUCTURED_GRID" << std::endl;

    output << "POINTS " << my_mesh->GetNnodes() << " float\n";
    for (int i = 0; i < my_mesh->GetNnodes(); i++)
    {
        auto node = std::dynamic_pointer_cast<ChNodeFEAxyz>(my_mesh->GetNode(i));
        output << node->GetPos().x() << " " << node->GetPos().y() << " " << node->GetPos().z() << "\n";
    }

    std::ifstream CopyFrom(MeshFileBuffer);
    output << CopyFrom.rdbuf();

    //Note that this is
    output << "\nPOINT_DATA " << my_mesh->GetNnodes() << "\n";
    output << "scalars strain float\n";
    output << "LOOKUP_TABLE default\n";
    for (int i = 0; i < my_mesh->GetNnodes(); i++)
    {
        double volAve = 0;
        ChVector<> StrainV;
        ChMatrixDynamic<> disp;
        double myVol = 0;
        double dx, dy, dz;
        for (int j = 0; j < NodeNeighborElement[i].size(); j++)
        {
            int myelemInx = NodeNeighborElement[i][j];
            // std::dynamic_pointer_cast<ChElementCableANCF>(my_mesh->GetElement(myelemInx))
            //     ->EvaluateSectionStrain(0.0, disp, StrainV);
            if (auto element = std::dynamic_pointer_cast<ChElementBrick>(my_mesh->GetElement(myelemInx)))
            {
                dx = std::dynamic_pointer_cast<ChElementBrick>(my_mesh->GetElement(myelemInx))->GetLengthX();
                dy = std::dynamic_pointer_cast<ChElementBrick>(my_mesh->GetElement(myelemInx))->GetLengthY();
                dz = std::dynamic_pointer_cast<ChElementBrick>(my_mesh->GetElement(myelemInx))->GetLengthZ();
                myVol += dx * dy * dz / NodeNeighborElement[i].size();
                volAve += 1 * dx * dy * dz / NodeNeighborElement[i].size();
            }
        }

        output << volAve / myVol << "\n";
    }


    output << "VECTORS Velocity float\n";
    for (int i = 0; i < my_mesh->GetNnodes(); i++) {
        ChVector<> vel = std::dynamic_pointer_cast<ChNodeFEAxyz>(my_mesh->GetNode(i))->GetPos_dt();
        vel += ChVector<>(1e-20);
        output << (double)vel.x() << " " << (double)vel.y() << " " << (double)vel.z() << "\n";
    }

    output << "VECTORS Acceleration float\n";

    for (int i = 0; i < my_mesh->GetNnodes(); i++) {
        ChVector<> acc = std::dynamic_pointer_cast<ChNodeFEAxyz>(my_mesh->GetNode(i))->GetPos_dtdt();
        acc += ChVector<>(1e-20);
        output << (double)acc.x() << " " << (double)acc.y() << " " << (double)acc.z() << "\n";
    }



    output.close();
};

#endif
