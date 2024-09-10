//
// Created by Giovane Avancini on 21/08/24.
//
#ifdef HAVE_CONFIG_H
#include <pz_config.h>
#endif

#include "meshpath_config.h"
#include <TPZGeoMeshTools.h>
#include <TPZAnalyticSolution.h>
#include <TPZGmshReader.h>
#include <TPZCompMeshTools.h>
#include <pzlog.h>
#include <pzshapequad.h>
#include <pzshapetriang.h>
#include <pzshapecube.h>
#include <pzshapetetra.h>
#include <TPZTimer.h>
#include <fstream>
#include <TPZSimpleTimer.h>
#include <TPZVTKGenerator.h>
#include <TPZRefPattern.h>
#include <tpzgeoelrefpattern.h>
#include <TPZRefPatternDataBase.h>
#include <pzbuildmultiphysicsmesh.h>
#include <pzstepsolver.h>       //for TPZStepSolver
#include <pzblockdiag.h>
#include <pzbdstrmatrix.h>
#include <pzcmesh.h>
#include <TPZMultiphysicsCompMesh.h>
#include <TPZLinearAnalysis.h>
#include <TPZVTKGeoMesh.h>
#include <pzskylstrmatrix.h>
#include <TPZYSMPMatrix.h>
#include <TPZNullMaterial.h>
#include <TPZNullMaterialCS.h>
#include <pzelementgroup.h>
#ifdef USING_MKL
#include <TPZYSMPPardiso.h>
#include <TPZSYSMPPardiso.h>
#include <TPZSSpStructMatrix.h> //symmetric sparse matrix storage
#endif
#include <pzcondensedcompel.h>
#include <pzelchdiv.h>
#include "ProblemData.h"
#include "TPZHybridCompressibleStokes.h"
#include "TPZInterfaceStokes.h"

const int global_nthread = 8;

TPZGeoMesh *CreateGMesh(ProblemData *inputData);

void InsertLagrangeMultipliers(ProblemData *inputData, TPZGeoMesh *gmesh);

void InsertInterfaces(TPZMultiphysicsCompMesh *cmesh, ProblemData *inputData, TPZGeoMesh *gmesh);

TPZCompMesh *CreateCMeshV(ProblemData *inputData, TPZGeoMesh *gmesh);

TPZCompMesh *CreateCMeshP(ProblemData *inputData, TPZGeoMesh *gmesh);

TPZCompMesh *CreateCMeshG(ProblemData *inputData, TPZGeoMesh *gmesh);

TPZCompMesh *CreateCMeshPm(ProblemData *inputData, TPZGeoMesh *gmesh);

TPZMultiphysicsCompMesh *CreateMultiPhysicsMesh(ProblemData *inputData, TPZGeoMesh *gmesh, TPZAnalyticSolution *sol, bool useIterative, REAL alpha = 0.1);

void CondenseElements(ProblemData *inputData, TPZMultiphysicsCompMesh *cmesh, TPZGeoMesh *gmesh, bool condensePressure = false);

void SolveProblemDirect(TPZLinearAnalysis &an, TPZCompMesh *cmesh, std::ofstream &outfile);

void SolveProblemIterative(TPZLinearAnalysis &an, TPZCompMesh *cmesh, REAL alpha, REAL tol, std::ofstream &outfile);

void PrintResults(TPZLinearAnalysis &an, TPZCompMesh *cmesh, int resolution = 0);

void IterativeSolver(TPZGeoMesh *gmesh, ProblemData* inputData, REAL alpha, TPZAnalyticSolution* sol, std::ofstream &outfile);

void DirectSolver(TPZGeoMesh *gmesh, ProblemData* inputData, TPZAnalyticSolution* sol, std::ofstream &outfile);

// Analytical solution
// constexpr int solOrder{4};
// auto exactSol = [](const TPZVec<REAL> &loc,
//                    TPZVec<STATE> &u,
//                    TPZFMatrix<STATE> &gradU)
// {
//     const auto &x = loc[0];
//     const auto &y = loc[1];
//     const auto &z = loc[2];

//     REAL aux = 1. / sinh(sqrt(2) * M_PI);
//     u[0] = sin(M_PI * x) * sin(M_PI * y) * sinh(sqrt(2) * M_PI * z) * aux;
//     gradU(0, 0) = M_PI * cos(M_PI * x) * sin(M_PI * y) * sinh(sqrt(2) * M_PI * z) * aux;
//     gradU(1, 0) = M_PI * cos(M_PI * y) * sin(M_PI * x) * sinh(sqrt(2) * M_PI * z) * aux;
//     gradU(2, 0) = sqrt(2) * M_PI * cosh(sqrt(2) * M_PI * z) * sin(M_PI * x) * sin(M_PI * y) * aux;
// };

int main(int argc, char *argv[])
{
#ifdef PZ_LOG
    TPZLogger::InitializePZLOG(std::string(MESHES_DIR) + "/" + "log4cxx.cfg");
#endif

    // Inputs
    ProblemData input;
    std::string filename = "Poiseuille";
    input.ReadJson(std::string(MESHES_DIR) + "/" + filename + ".json");

    const int ndiv = argc > 1 ? atoi(argv[1]) : 1;
    input.SetMeshName(std::string(MESHES_DIR) + "/" + filename + std::to_string(ndiv) + ".msh");
    if (argc > 2) //set p order at runtime
    {
        input.SetVelpOrder(atoi(argv[2])); 
        input.SetTracpOrder(atoi(argv[2])-1);
    }
    bool useIterative = argc > 3 ? atoi(argv[3]) : 1;
    REAL alpha = argc > 4 ? atof(argv[4]) : 0.001;

    // analytical solution
    TStokesAnalytic *exactSol = new TStokesAnalytic();
    exactSol->fvisco = input.DomainVec()[0].viscosity;
    exactSol->fvelocity = 1.;
    exactSol->fconstPressure = -0.5;
    exactSol->fDimension = input.Dim();
    exactSol->fExactSol = TStokesAnalytic::EPoisFlow;

    // geometric mesh
    TPZGeoMesh *gmesh = CreateGMesh(&input);
    InsertLagrangeMultipliers(&input, gmesh);
    {
        std::ofstream out("gmesh.vtk");
        TPZVTKGeoMesh::PrintGMeshVTK(gmesh, out);
    }

    std::string output_name = "stokes-ndiv-" + std::to_string(ndiv) + "-p-" + std::to_string(input.VelpOrder()) + "-iter-" + std::to_string(useIterative);
    if (argc > 4)
        output_name += "-alpha-" + std::string(argv[4]);
    output_name += ".dat";
    std::ofstream outfile(output_name);

    if (useIterative)
    {
        IterativeSolver(gmesh, &input, alpha, exactSol, outfile);
    }
    else
    {
        DirectSolver(gmesh, &input, exactSol, outfile);
    }

    return 0;
}

TPZGeoMesh *CreateGMesh(ProblemData *inputData)
{
    TPZGeoMesh *gmesh = new TPZGeoMesh;
    gmesh->SetName("GeoMesh");
    TPZGmshReader reader;
    reader.GeometricGmshMesh(inputData->MeshName(), gmesh);
    gmesh->BuildConnectivity();

    return gmesh;
}

void InsertLagrangeMultipliers(ProblemData *inputData, TPZGeoMesh *gmesh)
{
    int64_t nEl = gmesh->NElements();

    // we look for tangential BCs
    TPZVec<int> IDVec(inputData->TangentialBCs().size(), 0);

    for (int i = 0; i < inputData->TangentialBCs().size(); i++)
        IDVec[i] = inputData->TangentialBCs()[i].matID;

    for (auto const &BcMatID : IDVec)
    {
        for (int64_t el = 0; el < nEl; el++)
        {
            int meshDim = gmesh->Dimension();

            TPZGeoEl *geoEl = gmesh->Element(el);
            int matID = geoEl->MaterialId();

            if (matID != BcMatID)
                continue;

            int nSides = geoEl->NSides();
            TPZGeoElSide geoElSide(geoEl, nSides - 1);
            TPZCompElSide compElSide = geoElSide.Reference();

            TPZStack<TPZGeoElSide> neighbourSet;
            geoElSide.AllNeighbours(neighbourSet);

            int64_t nneighs = neighbourSet.size();

            for (int stack_i = 0; stack_i < nneighs; stack_i++)
            {
                TPZGeoElSide neighbour = neighbourSet[stack_i];
                int neighMatID = neighbour.Element()->MaterialId();
                TPZCompElSide compElNeigh = neighbour.Reference();

                int64_t neighIndex = neighbour.Element()->Index();

                if (neighbour.Element()->Dimension() != meshDim)
                    continue;

                if (neighbour.Element()->HasSubElement())
                    DebugStop();

                TPZGeoElBC(neighbour, inputData->InterfaceID());

                neighbour = neighbour.Neighbour();

                if (neighbour.Element()->MaterialId() != inputData->InterfaceID())
                    DebugStop();

                TPZGeoElBC(neighbour, inputData->LambdaID());

                neighbour = neighbour.Neighbour();

                if (neighbour.Element()->MaterialId() != inputData->LambdaID())
                    DebugStop();

                TPZGeoElBC(neighbour, 21);
            }
        }
    }

    // we look for two domain neighbour elements
    for (int64_t el = 0; el < nEl; el++)
    {
        TPZGeoEl *geoEl = gmesh->Element(el);

        if (!geoEl)
            continue;
        if (geoEl->HasSubElement())
            continue;
        if (geoEl->Dimension() != gmesh->Dimension())
            continue;

        int nside = geoEl->NSides();

        if (inputData->DomainVec().size() != 0)
        {
            for (int side = 0; side < nside; side++)
            {
                if (geoEl->SideDimension(side) != gmesh->Dimension() - 1)
                    continue;

                TPZGeoElSide geoElSide(geoEl, side);
                TPZGeoElSide neighbour = geoElSide.Neighbour();

                if (neighbour == geoElSide)
                    continue;
                if (neighbour.Element()->HasSubElement())
                    continue;

                TPZGeoElSide neighbour2 = neighbour;
                while (neighbour2 != geoElSide)
                {
                    if (neighbour2.Element()->MaterialId() == inputData->InterfaceID())
                    {
                        break;
                    }

                    if (neighbour2.Element()->Dimension() == gmesh->Dimension())
                        neighbour = neighbour2;

                    neighbour2 = neighbour2.Neighbour();
                }

                if (neighbour2 == geoElSide)
                {
                    TPZGeoElBC(geoElSide, inputData->InterfaceID());
                    TPZGeoElBC(neighbour, inputData->InterfaceID());

                    neighbour2 = neighbour.Neighbour();
                    neighbour = geoElSide.Neighbour();

                    if (neighbour.Element()->MaterialId() != inputData->InterfaceID() || neighbour2.Element()->MaterialId() != inputData->InterfaceID())
                        DebugStop();

                    TPZGeoElBC(neighbour, inputData->LambdaID());
                    TPZGeoElBC(neighbour2, inputData->LambdaID());

                    neighbour = neighbour.Neighbour();
                    neighbour2 = neighbour2.Neighbour();

                    if (neighbour.Element()->MaterialId() != inputData->LambdaID() || neighbour2.Element()->MaterialId() != inputData->LambdaID())
                        DebugStop();

                    TPZGeoElBC(neighbour, 21);
                    TPZGeoElBC(neighbour2, 21);

                    neighbour = neighbour.Neighbour();

                    if (neighbour.Element()->MaterialId() != 21)
                        DebugStop();

                    TPZGeoElBC(neighbour, inputData->TanVelID());
                }
            }
        }
    }

    gmesh->BuildConnectivity();
}

void InsertInterfaces(TPZMultiphysicsCompMesh *cmesh, ProblemData *inputData, TPZGeoMesh *gmesh)
{
    if (!gmesh)
        DebugStop();

    TPZManVector<int64_t, 3> LeftElIndices(1, 0), RightElIndices(1, 1);

    int nInterfaceCreated = 0;

    int mat_lambda = inputData->LambdaID();
    int mat_tan = inputData->TanVelID();

    int meshDim = cmesh->Dimension();

    std::set<int> BcIDs;
    for (int i = 0; i < inputData->TangentialBCs().size(); i++)
        BcIDs.insert(inputData->TangentialBCs()[i].matID);

    int64_t nEl = gmesh->NElements();

    for (int64_t el = 0; el < nEl; el++)
    {
        TPZGeoEl *geoEl = gmesh->Element(el);
        int matID = geoEl->MaterialId();

        if (geoEl->Dimension() != meshDim)
            continue;

        if (geoEl->HasSubElement())
            continue;

        int nSides = geoEl->NSides();

        for (int side = 0; side < nSides; side++)
        {
            if (geoEl->SideDimension(side) != meshDim - 1)
                continue;

            TPZGeoElSide geoElSide(geoEl, side);
            TPZGeoElSide geoEl_neighbour = geoElSide.Neighbour();
            int neighbour_matID = geoEl_neighbour.Element()->MaterialId();

            if (neighbour_matID != inputData->InterfaceID())
            {
                int a = 0;
            }

            if (neighbour_matID == inputData->InterfaceID())
            {
                TPZCompElSide compElSide = geoElSide.Reference();

                if (!compElSide)
                    DebugStop();

                TPZGeoElSide geoEl_neighbour2 = geoEl_neighbour.Neighbour();

                if (geoEl_neighbour2.Element()->MaterialId() == inputData->ObstructionID())
                    DebugStop();

                if (geoEl_neighbour2.Element()->MaterialId() != inputData->LambdaID())
                    continue;

                TPZCompElSide compEl_neighbour = geoEl_neighbour2.Reference();

                if (!compEl_neighbour)
                    DebugStop();

                if (geoEl_neighbour2.Element()->HasSubElement())
                    DebugStop();

                // creating the interface between the Hdiv domain element and tangential stress
                TPZMultiphysicsInterfaceElement *interElem = new TPZMultiphysicsInterfaceElement(*cmesh, geoEl_neighbour.Element(), compElSide, compEl_neighbour);
                interElem->SetLeftRightElementIndices(LeftElIndices, RightElIndices);
                nInterfaceCreated++;

                // we take the second interface between tangential stress and tangential velocity
                geoEl_neighbour = geoEl_neighbour2.Neighbour();

                if (geoEl_neighbour.Element()->MaterialId() != 21)
                    DebugStop();

                TPZGeoElSide geoEl_neighbour3 = geoEl_neighbour;
                for (; geoEl_neighbour3 != geoElSide; geoEl_neighbour3++)
                {
                    int neighbour_matID2 = geoEl_neighbour3.Element()->MaterialId();

                    if (neighbour_matID2 != mat_tan && BcIDs.find(neighbour_matID2) == BcIDs.end())
                        continue;

                    if (neighbour_matID2 == inputData->ObstructionID())
                        DebugStop();

                    compElSide = geoEl_neighbour3.Reference();

                    if (!compElSide)
                        DebugStop();

                    if (geoEl_neighbour3.Element()->HasSubElement())
                        DebugStop();

                    TPZMultiphysicsInterfaceElement *interElem = new TPZMultiphysicsInterfaceElement(*cmesh, geoEl_neighbour.Element(), compElSide, compEl_neighbour);
                    interElem->SetLeftRightElementIndices(LeftElIndices, RightElIndices);
                    nInterfaceCreated++;

                    break;
                }
            }
        }
    }

    std::cout << __PRETTY_FUNCTION__ << "Number of Interfaces Created: " << nInterfaceCreated << std::endl;
}

TPZCompMesh *CreateCMeshV(ProblemData *inputData, TPZGeoMesh *gmesh)
{
    TPZCompMesh *cmesh_v = new TPZCompMesh(gmesh);
    cmesh_v->SetName("Cmesh_V");

    std::set<int> materialID;

    cmesh_v->ApproxSpace().CreateDisconnectedElements(false);

    // domain's material (2D or 3D)
    if (inputData->DomainVec().size() != 0)
    {
        cmesh_v->SetDefaultOrder(inputData->VelpOrder());
        cmesh_v->SetDimModel(inputData->Dim());

        if (inputData->HdivType() == ProblemData::EConstant)
        {
            cmesh_v->ApproxSpace().SetHDivFamily(HDivFamily::EHDivConstant);
        }
        else if (inputData->HdivType() == ProblemData::EStandard)
        {
            cmesh_v->ApproxSpace().SetHDivFamily(HDivFamily::EHDivStandard);
        }

        cmesh_v->SetAllCreateFunctionsHDiv();

        auto *mat_normal = new TPZNullMaterial<>(inputData->DomainVec()[0].matID);
        cmesh_v->InsertMaterialObject(mat_normal);

        materialID.insert(inputData->DomainVec()[0].matID);

        // boundary condition's material
        TPZFMatrix<STATE> val1(1, 1, 0.);
        TPZManVector<STATE> val2(1, 0.);

        for (const auto &bc : inputData->NormalBCs())
        {
            val2 = bc.value;

            auto BCmat = mat_normal->CreateBC(mat_normal, bc.matID, bc.type, val1, val2);
            cmesh_v->InsertMaterialObject(BCmat);
            materialID.insert(bc.matID);
        }

        cmesh_v->AutoBuild(materialID);

        // Increasing internal function order
        int64_t ncEl = cmesh_v->NElements();

        for (int64_t cEl = 0; cEl < ncEl; cEl++)
        {
            TPZCompEl *compEl = cmesh_v->Element(cEl);

            if (compEl->Dimension() == inputData->Dim())
            {
                TPZInterpolatedElement *intercEl = dynamic_cast<TPZInterpolatedElement *>(compEl);

                // checking of the dynamic cast existis
                if (!intercEl)
                    continue;

                intercEl->ForceSideOrder(compEl->Reference()->NSides() - 1, inputData->VelpOrder() + 2);
            }
        }

        int64_t nconV = cmesh_v->NConnects();

        gmesh->ResetReference();
        materialID.clear();

        // tangential velocity material
        auto mat_tan = new TPZNullMaterial<>(inputData->TanVelID());
        mat_tan->SetNStateVariables(inputData->Dim() - 1); // in 3d, there are 2 state variables (one at each tangential direction)
        cmesh_v->InsertMaterialObject(mat_tan);

        materialID.insert(inputData->TanVelID());

        // tangential velocity on boundary material
        for (const auto &bc : inputData->TangentialBCs())
        {
            auto matBC = new TPZNullMaterial<>(bc.matID);
            matBC->SetNStateVariables(inputData->Dim() - 1);
            cmesh_v->InsertMaterialObject(matBC);

            materialID.insert(bc.matID);
        }

        cmesh_v->ApproxSpace().CreateDisconnectedElements(true);
        cmesh_v->SetDefaultOrder(inputData->TracpOrder());
        cmesh_v->SetDimModel(inputData->Dim() - 1);
        cmesh_v->AutoBuild(materialID);

        int64_t ncon = cmesh_v->NConnects();
        for (int64_t i = nconV; i < ncon; i++)
        {
            TPZConnect &newnod = cmesh_v->ConnectVec()[i];
                newnod.SetLagrangeMultiplier(1);
        }

        gmesh->ResetReference();
    }

    cmesh_v->ExpandSolution();
    inputData->MeshVector()[0] = cmesh_v;

    return cmesh_v;
}

TPZCompMesh *CreateCMeshP(ProblemData *inputData, TPZGeoMesh *gmesh)
{
    TPZCompMesh *cmesh_p = new TPZCompMesh(gmesh);
    cmesh_p->SetName("CMesh_P");

    std::set<int> materialID;

    if (inputData->DomainVec().size() != 0)
    {
        cmesh_p->SetDimModel(inputData->Dim());

        if (inputData->HdivType() == ProblemData::EConstant)
        {
            cmesh_p->SetAllCreateFunctionsDiscontinuous();
            cmesh_p->SetDefaultOrder(0);
        }
        else if (inputData->HdivType() == ProblemData::EStandard)
        {
            cmesh_p->SetDefaultOrder(inputData->VelpOrder() + 2);
            cmesh_p->SetAllCreateFunctionsContinuous();
        }

        cmesh_p->ApproxSpace().CreateDisconnectedElements(true);

        // domain's material
        auto *mat = new TPZNullMaterial<>(inputData->DomainVec()[0].matID);
        cmesh_p->InsertMaterialObject(mat);

        materialID.insert(inputData->DomainVec()[0].matID);

        cmesh_p->AutoBuild(materialID);
        gmesh->ResetReference();

        materialID.clear();

        int64_t ncon = cmesh_p->NConnects();
        for (int64_t i = 0; i < ncon; i++)
        {
            TPZConnect &newnod = cmesh_p->ConnectVec()[i];
            newnod.SetLagrangeMultiplier(2);
        }

        // matLambda traction material
        auto matLambda = new TPZNullMaterial<>(inputData->LambdaID());
        matLambda->SetNStateVariables(inputData->Dim() - 1);
        cmesh_p->InsertMaterialObject(matLambda);

        materialID.insert(inputData->LambdaID());

        if (inputData->TracpOrder() > 0)
        {
            cmesh_p->SetAllCreateFunctionsContinuous();
            cmesh_p->ApproxSpace().CreateDisconnectedElements(true);
        }
        else
        {
            cmesh_p->SetAllCreateFunctionsDiscontinuous();
            cmesh_p->ApproxSpace().CreateDisconnectedElements(true);
        }

        cmesh_p->SetDefaultOrder(inputData->TracpOrder());
        cmesh_p->SetDimModel(inputData->Dim() - 1);
        cmesh_p->AutoBuild(materialID);

        gmesh->ResetReference();

        ncon = cmesh_p->NConnects();
        for (int64_t i = 0; i < ncon; i++)
        {
            TPZConnect &newnod = cmesh_p->ConnectVec()[i];

            if (newnod.LagrangeMultiplier() == 0)
                newnod.SetLagrangeMultiplier(3);
        }
    }

    cmesh_p->ExpandSolution();
    inputData->MeshVector()[1] = cmesh_p;

    return cmesh_p;
}

TPZCompMesh *CreateCMeshG(ProblemData *inputData, TPZGeoMesh *gmesh)
{
    TPZCompMesh *cmesh_g = new TPZCompMesh(gmesh);

    cmesh_g->SetName("CMesh_g");
    cmesh_g->SetDefaultOrder(0);
    cmesh_g->SetDimModel(inputData->Dim());

    cmesh_g->SetAllCreateFunctionsDiscontinuous();

    auto mat = new TPZNullMaterial<>(inputData->DomainVec()[0].matID);
    cmesh_g->InsertMaterialObject(mat);

    cmesh_g->AutoBuild();
    cmesh_g->AdjustBoundaryElements();
    cmesh_g->CleanUpUnconnectedNodes();

    int64_t ncon = cmesh_g->NConnects();
    for (int64_t i = 0; i < ncon; i++)
    {
        TPZConnect &newnod = cmesh_g->ConnectVec()[i];
        newnod.SetLagrangeMultiplier(4);
    }

    inputData->MeshVector()[2] = cmesh_g;

    return cmesh_g;
}

TPZCompMesh *CreateCMeshPm(ProblemData *inputData, TPZGeoMesh *gmesh)
{
    TPZCompMesh *cmesh_pm = new TPZCompMesh(gmesh);

    cmesh_pm->SetName("CMesh_pm");
    cmesh_pm->SetDefaultOrder(0);
    cmesh_pm->SetDimModel(inputData->Dim());

    cmesh_pm->SetAllCreateFunctionsDiscontinuous();

    auto mat = new TPZNullMaterial<>(inputData->DomainVec()[0].matID);
    cmesh_pm->InsertMaterialObject(mat);

    cmesh_pm->AutoBuild();
    cmesh_pm->AdjustBoundaryElements();
    cmesh_pm->CleanUpUnconnectedNodes();

    int64_t ncon = cmesh_pm->NElements();

    for (int64_t i = 0; i < ncon; i++)
    {
        TPZConnect &newnod = cmesh_pm->ConnectVec()[i];
        newnod.SetLagrangeMultiplier(5);
    }

    inputData->MeshVector()[3] = cmesh_pm;

    return cmesh_pm;
}

TPZMultiphysicsCompMesh *CreateMultiPhysicsMesh(ProblemData *inputData, TPZGeoMesh *gmesh, TPZAnalyticSolution *sol, bool useIterative, REAL alpha)
{
    TPZMultiphysicsCompMesh *cmesh_m = new TPZMultiphysicsCompMesh(gmesh);

    cmesh_m->SetName("CMesh_M");

    cmesh_m->SetDefaultOrder(inputData->VelpOrder());
    cmesh_m->SetAllCreateFunctionsMultiphysicElem();

    // Creating Materials
    std::set<int> materialID;

    // 1. For domain
    if (inputData->DomainVec().size() != 0)
    {
        REAL viscosity = inputData->DomainVec()[0].viscosity;

        if (dynamic_cast<TStokesAnalytic *>(sol))
        {
            TStokesAnalytic *flow = dynamic_cast<TStokesAnalytic *>(sol);
            viscosity = flow->fvisco;

            if (flow->fExactSol == TStokesAnalytic::ENone)
                sol = nullptr;
        }

        TPZHybridStokes *domain_mat;
        if (useIterative)
            domain_mat = new TPZHybridCompressibleStokes(inputData->DomainVec()[0].matID, inputData->Dim(), viscosity, alpha);
        else
            domain_mat = new TPZHybridStokes(inputData->DomainVec()[0].matID, inputData->Dim(), viscosity);
        if (sol)
        {
            domain_mat->SetExactSol(sol->ExactSolution(), 3);
            domain_mat->SetForcingFunction(sol->ForceFunc(), 3);
        }

        cmesh_m->InsertMaterialObject(domain_mat);
        materialID.insert(inputData->DomainVec()[0].matID);

        // 2. Boundary Conditions
        TPZFMatrix<STATE> val1(3, 3, 0.);
        TPZManVector<STATE> val2(3, 0.);

        // Normal
        for (const auto &bc : inputData->NormalBCs())
        {
            val2 = bc.value;

            TPZBndCond *matBC = domain_mat->CreateBC(domain_mat, bc.matID, bc.type, val1, val2);

            auto matBC2 = dynamic_cast<TPZBndCondT<STATE> *>(matBC);
            if (sol)
                matBC2->SetForcingFunctionBC(sol->ExactSolution(), 3);

            cmesh_m->InsertMaterialObject(matBC);
            materialID.insert(bc.matID);
        }

        // Tangential
        for (const auto &bc : inputData->TangentialBCs())
        {
            val2 = bc.value;

            TPZBndCond *matBC = domain_mat->CreateBC(domain_mat, bc.matID, bc.type, val1, val2);

            auto matBC2 = dynamic_cast<TPZBndCondT<STATE> *>(matBC);
            if (sol)
                matBC2->SetForcingFunctionBC(sol->ExactSolution(), 3);

            cmesh_m->InsertMaterialObject(matBC);
            materialID.insert(bc.matID);
        }

        // 3. Material for Tangential Velocity
        TPZNullMaterialCS<> *mat_tan = new TPZNullMaterialCS<>(inputData->TanVelID());
        mat_tan->SetDimension(inputData->Dim() - 1);
        mat_tan->SetNStateVariables(inputData->Dim() - 1);

        cmesh_m->InsertMaterialObject(mat_tan);
        materialID.insert(inputData->TanVelID());

        // 4. Material for Tangential Traction
        TPZNullMaterialCS<> *mat_lambda = new TPZNullMaterialCS<>(inputData->LambdaID());
        mat_lambda->SetDimension(inputData->Dim() - 1);
        mat_lambda->SetNStateVariables(inputData->Dim() - 1);

        cmesh_m->InsertMaterialObject(mat_lambda);
        materialID.insert(inputData->LambdaID());
    }

    TPZManVector<int, 2> active_approx_spaces(inputData->MeshVector().size(), 1);

    cmesh_m->BuildMultiphysicsSpace(active_approx_spaces, inputData->MeshVector());
    cmesh_m->AdjustBoundaryElements();
    cmesh_m->CleanUpUnconnectedNodes();
    cmesh_m->LoadReferences();

    // 5. Material for Interfaces (inner)
    TPZInterfaceStokes *matInterface = new TPZInterfaceStokes(inputData->InterfaceID(), inputData->Dim() - 1);
    matInterface->SetMultiplier(1.);

    cmesh_m->InsertMaterialObject(matInterface);

    TPZInterfaceStokes *matInterface2 = new TPZInterfaceStokes(21, inputData->Dim() - 1);
    matInterface2->SetMultiplier(1.);

    cmesh_m->InsertMaterialObject(matInterface2);

    InsertInterfaces(cmesh_m, inputData, gmesh);

    if (inputData->CondensedElements())
    {
        cmesh_m->SetName("CMesh_M_BeforeCond");
        cmesh_m->ComputeNodElCon();
    }

    return cmesh_m;
}

void CondenseElements(ProblemData *inputData, TPZMultiphysicsCompMesh *cmesh_m, TPZGeoMesh *gmesh, bool condensePressure)
{
    int64_t nCompEl = cmesh_m->ElementVec().NElements();
    int dim = gmesh->Dimension();

    std::set<int64_t> externalNodes;
    TPZStack<int64_t> groupIndex;
    TPZStack<TPZElementGroup*> elGroups;
    int count = 0;

    std::set<int> BCsIDs;
    for (int i = 0; i < inputData->TangentialBCs().size(); i++)
        BCsIDs.insert(inputData->TangentialBCs()[i].matID);
    for (int i = 0; i < inputData->NormalBCs().size(); i++)
        BCsIDs.insert(inputData->NormalBCs()[i].matID);

    // Creating the element groups for the domain
    for (int64_t el = 0; el < nCompEl; el++)
    {
        TPZCompEl *compEl = cmesh_m->Element(el);

        if (compEl->Dimension() != dim)
            continue;

        int nConnect = compEl->NConnects();

        if (!condensePressure)
        {
            int64_t coIndex = compEl->ConnectIndex(nConnect - 1);
            externalNodes.insert(coIndex);
        }

        count++;
        groupIndex.push_back(compEl->Index());

        TPZElementGroup *groupEl = new TPZElementGroup(*cmesh_m);
        elGroups.Push(groupEl);
        elGroups[count - 1]->AddElement(compEl);
    }

    int64_t nGeoEl = gmesh->NElements();
    for (int elgr = 0; elgr < elGroups.size(); elgr++)
    {
        TPZElementGroup *groupEl = elGroups[elgr];
        TPZCompEl *compEl = groupEl->GetElGroup()[0];

        if (!compEl)
            DebugStop();

        TPZGeoEl *geoEl = compEl->Reference();

        if (geoEl->Dimension() != dim)
            DebugStop();

        if (geoEl->HasSubElement())
            DebugStop();

        int nSides = geoEl->NSides();

        int64_t compEl_Index = compEl->Index();

        for (int side = 0; side < nSides; side++)
        {
            if (geoEl->SideDimension(side) != dim - 1)
                continue;

            TPZGeoElSide geoEl_Side(geoEl, side);

            TPZStack<TPZGeoElSide> allNeighbours;
            geoEl_Side.AllNeighbours(allNeighbours);

            TPZStack<TPZCompEl *> neighbourCompEls;

            if (allNeighbours.size() > 0)
            {

                if (allNeighbours[0].Element()->MaterialId() == inputData->InterfaceID())
                {
                    for (int i = 0; i < 3; i++)
                    {
                        TPZGeoEl *gel = allNeighbours[i].Element();
                        int matID = gel->MaterialId();

                        if (matID != inputData->InterfaceID() && matID != inputData->LambdaID() && matID != 21)
                            DebugStop();

                        TPZCompEl *cel = allNeighbours[i].Element()->Reference();

                        if (!cel)
                            DebugStop();

                        groupEl->AddElement(cel);
                    }
                }
            }
        }
    }

    cmesh_m->ComputeNodElCon();

    for (auto it = externalNodes.begin(); it != externalNodes.end(); it++)
    {
        int64_t coIndex = *it;
        cmesh_m->ConnectVec()[coIndex].IncrementElConnected();
    }

    // Creating condensed elements
    int64_t nEnvel = elGroups.NElements();
    for (int64_t iEnv = 0; iEnv < nEnvel; iEnv++)
    {
        TPZElementGroup *elGroup = elGroups[iEnv];
        new TPZCondensedCompElT<STATE>(elGroup, false);
    }

    cmesh_m->SetName("CMesh_M_Condensed");
    cmesh_m->CleanUpUnconnectedNodes();
    cmesh_m->ExpandSolution();
}

void SolveProblemDirect(TPZLinearAnalysis &an, TPZCompMesh *cmesh, std::ofstream &outfile)
{
    #ifdef USING_MKL
    TPZSSpStructMatrix<STATE,TPZStructMatrixOR<STATE>> strmat(cmesh);
#else
    TPZSkylineStructMatrix<STATE> strmat(cmesh);
#endif
    strmat.SetNumThreads(global_nthread);

    an.SetStructuralMatrix(strmat);

    /// Setting a direct solver
    TPZStepSolver<STATE> step;
    step.SetDirect(ELDLt); // ELU //ECholesky // ELDLt

    an.SetSolver(step);

    // assembles the system
    std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
    an.Assemble();
    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
    std::cout << "Time Assemble = " << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << "[ms]" << std::endl;
    if (!outfile.fail())
        outfile << "Time Assemble = " << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << "[ms]" << std::endl;

    /// solves the system
    std::chrono::steady_clock::time_point begin2 = std::chrono::steady_clock::now();
    an.Solve();
    std::chrono::steady_clock::time_point end2 = std::chrono::steady_clock::now();
    std::cout << "Time Solve = " << std::chrono::duration_cast<std::chrono::milliseconds>(end2 - begin2).count() << "[ms]" << std::endl;
    if (!outfile.fail())
        outfile << "Time Solve Direct = " << std::chrono::duration_cast<std::chrono::milliseconds>(end2 - begin2).count() << "[ms]" << std::endl;
}

void SolveProblemIterative(TPZLinearAnalysis &an, TPZCompMesh *cmesh, REAL alpha, REAL tol, std::ofstream &outfile)
{
#ifdef USING_MKL
    TPZSSpStructMatrix<STATE,TPZStructMatrixOR<STATE>> strmat(cmesh);
#else
    TPZSkylineStructMatrix<STATE> strmat(cmesh);
#endif
    strmat.SetNumThreads(global_nthread);
    an.SetStructuralMatrix(strmat);

    /// Setting a direct solver
    TPZStepSolver<STATE> step;
    step.SetDirect(ECholesky); // ELU //ECholesky // ELDLt

    an.SetSolver(step);

    // assembles the system
    std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
    an.Assemble();
    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
    std::cout << "Time Assemble = " << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << "[ms]" << std::endl;
    if (!outfile.fail())
        outfile << "Time Assemble = " << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << "[ms]" << std::endl;

    // Computing the number of hdiv domain elements, the number of facet connects and the number of facet equations
    TPZMultiphysicsCompMesh *cmesh_m = dynamic_cast<TPZMultiphysicsCompMesh *>(cmesh);
    TPZCompMesh *cmesh_u = cmesh_m->MeshVector()[0];
    int64_t nel_u = cmesh_u->NElements();
    int64_t nElements = 0;
    int64_t nFacetConnects = 0;
    int64_t nFacetEqs = 0;

    for (int64_t iel = 0; iel < nel_u; iel++)
    {
        TPZCompEl *cel = cmesh_u->Element(iel);
        if (!cel)
            continue;
        if (cel->Dimension() != cmesh->Dimension())
            continue;

        TPZGeoEl *gel = cel->Reference();
        nFacetConnects += gel->NSides(cmesh->Dimension() - 1);
        nElements++;
    }

    for (auto &c : cmesh_u->ConnectVec())
    {
        if (c.LagrangeMultiplier() != 0 || c.HasDependency() || c.IsCondensed() || c.NElConnected() == 1)
            continue;
        nFacetEqs += c.NShape();
    }

    /*
    For Hdiv constant space, the global matrix B has dimensions (nFacetEqs x nElements)
    and only has contribution of the constant flux dof of each face equal to -fsideorient.
    Thus the number of nonzeros of B is nFaceConnects.
    */
    TPZVec<int64_t> iBT(nElements + 1, 0);
    TPZVec<int64_t> jBT(nFacetConnects, 0);
    TPZVec<STATE> valBT(nFacetConnects, 0.);
    TPZVec<STATE> elArea(nElements, 0.);

    int64_t row = 0;
    int count = 0;
    for (int64_t iel = 0; iel < nel_u; iel++)
    {
        int64_t col = 0;
        TPZCompEl *cel = cmesh_u->Element(iel);
        if (!cel)
            continue;
        if (cel->Dimension() != cmesh->Dimension())
            continue;
        TPZGeoEl *gel = cel->Reference();
        int nsides = gel->NSides();
        int nfacets = gel->NSides(gel->Dimension() - 1);
        int eqcont = 0;
        iBT[row + 1] = iBT[row] + nfacets;
        std::set<int64_t> connectlist;
        cel->BuildConnectList(connectlist);
        for (auto it = connectlist.begin(); it != connectlist.end(); it++)
        {
            TPZConnect &con = cmesh_m->ConnectVec()[*it];

            if (con.LagrangeMultiplier() != 0 || con.HasDependency() || con.IsCondensed() || con.NElConnected() == 1)
                continue;

            for (int iside = gel->NCornerNodes(); iside < nsides - 1; iside++)
            {
                if (gel->SideDimension(iside) != cmesh->Dimension() - 1)
                    continue; // only the facets
                TPZCompElSide celside(cel, iside);
                int64_t connect_id = celside.ConnectIndex();
                if (connect_id == *it)
                {
                    REAL side_orient = gel->NormalOrientation(iside);
                    valBT[count] = -side_orient;
                    break;
                }
            }
            int64_t seq = con.SequenceNumber();
            col = cmesh_m->Block().Position(seq);
            jBT[count] = col;
            count++;
        }
        const TPZIntPoints &intpoints = cel->GetIntegrationRule();
        int geldim = gel->Dimension();
        TPZManVector<REAL, 3> xp(3, 0.);
        int np = intpoints.NPoints();
        for (int ip = 0; ip < np; ip++)
        {
            TPZManVector<REAL, 3> xi(geldim, 0.);
            REAL weight;
            intpoints.Point(ip, xi, weight);
            gel->X(xi, xp);
            REAL detjac;
            TPZFMatrix<REAL> jac(geldim, geldim), jacinv(geldim, geldim);
            TPZFMatrix<REAL> axes(geldim, 3);
            REAL jacdet;
            gel->Jacobian(xp, jac, axes, jacdet, jacinv);
            elArea[row] += weight * jacdet;
        }
        row++;
    }

    // Gettig the Global stiffness matrix from solver
    auto KG = an.MatrixSolver<STATE>().Matrix();
    KG->SetDefPositive(true);
    int64_t neq = KG->Rows();
    TPZFYsmpMatrix<STATE> BT(nElements, neq);
    BT.SetData(iBT, jBT, valBT);


    TPZFMatrix<STATE> force = an.Rhs();
    TPZFMatrix<STATE> rhs(cmesh->NEquations(), 1, 0.);
    TPZFMatrix<STATE> dsol(cmesh->NEquations(), 1, 0.);

    // Obtaining the initial solution
    begin = std::chrono::steady_clock::now();
    an.Solve();
    TPZFMatrix<STATE> sol = an.Solution();
    end = std::chrono::steady_clock::now();
    std::cout << "Time for Initial Solver = " << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << "[ms]" << std::endl;
    if (!outfile.fail())
        outfile << "Time for Initial Solver = " << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << "[ms]" << std::endl;
    REAL timesolve = std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count();

    // Computing the residual (force is not used)
    TPZFMatrix<STATE> aux(nElements, 1, 0.);
    BT.MultAdd(sol, force, aux, 1.0, 0.0); // aux = BT * sol
    for (int64_t i = 0; i < nElements; i++)
    {
        aux(i, 0) /= elArea[i];
    }
    BT.MultAdd(aux, force, rhs, -1.0 / alpha, 0.0, 1); // rhs = -Transpose(BT) * aux / alpha = -Transpose(BT) * BT * sol / alpha
    dsol = rhs;

    // Computing initial pressure (force is not used)
    TPZFMatrix<STATE> pressure(nElements, 1, 0.);
    BT.MultAdd(sol, force, pressure, 1.0 / alpha, 0.0); // pressure = BT * sol / alpha
    for (int64_t i = 0; i < nElements; i++)
    {
        pressure(i, 0) /= elArea[i];
    }

    REAL norm_dsol = 1.0, norm_rhs = 1.0;
    int nit = 0;
    const int size = rhs.Rows();
    while (norm_dsol > tol || norm_rhs > tol)
    {
        begin = std::chrono::steady_clock::now();
        KG->SolveDirect(dsol, ECholesky);
        end = std::chrono::steady_clock::now();
        REAL iterativetime = std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count();
        timesolve += iterativetime;
        sol += dsol;
        // Updating the residual (force is not used)
        BT.MultAdd(sol, force, aux, 1.0, 0.0); // aux = BT * sol
        for (int64_t i = 0; i < nElements; i++)
        {
            aux(i, 0) /= elArea[i];
        }

        BT.MultAdd(aux, force, rhs, -1.0 / alpha, 0.0, 1); // rhs = Transpose(BT) * aux = Transpose(BT) * BT * sol
        // Updating pressure (force is not used)
        TPZFMatrix<STATE> dp(nElements, 1, 0.);
        BT.MultAdd(sol, force, dp, 1.0 / alpha, 0.0); // dp = BT * sol / alpha
        for (int64_t i = 0; i < nElements; i++)
        {
            dp(i, 0) /= elArea[i];
        }

        pressure += dp;
        nit++;
        norm_dsol = 0.0;
        norm_rhs = 0.0;
        for (int64_t i = 0; i < size; i++)
        {
            norm_dsol += dsol(i, 0) * dsol(i, 0);
            norm_rhs += rhs(i, 0) * rhs(i, 0);
        }
        norm_dsol = sqrt(norm_dsol);
        norm_rhs = sqrt(norm_rhs);
        dsol = rhs;
        std::cout << "Iteration: " << nit << ". Time spent: " << iterativetime << ". dsol_norm: " << norm_dsol << ", rhs_norm: " << norm_rhs << std::endl;
        if (!outfile.fail())
            outfile << "Iteration: " << nit << ". Time spent: " << iterativetime << ". dsol_norm: " << norm_dsol << ", rhs_norm: " << norm_rhs << std::endl;

        if (nit > 50)
        {
            std::cout << "Solver diverged.\n";
            break;
        }
    }
    std::cout << "Time iterative solver = " << timesolve << "[ms]" << std::endl;
    if (!outfile.fail())
        outfile << "Time iterative solver = " << timesolve << "[ms]" << std::endl;

    // Transfering the solution to the mesh
    an.Solution() = sol;
    an.LoadSolution();

    TPZFMatrix<REAL> &mesh_sol = cmesh->Solution();
    const int64_t nconnectHDiv = cmesh_u->NConnects();
    const int64_t nEquationsFull = mesh_sol.Rows();
    for (int64_t ic = nconnectHDiv; ic < nconnectHDiv+nElements; ic++)
    {
        TPZConnect &con = cmesh->ConnectVec()[ic];
        int64_t seq = con.SequenceNumber();
        int64_t pos = cmesh->Block().Position(seq);
        mesh_sol(pos, 0) = pressure(ic-nconnectHDiv, 0);
    }
    // for (const TPZConnect &con : cmesh->ConnectVec())
    // {
    //     if (con.LagrangeMultiplier() != 2)
    //         continue; // only pressure connect
    //     int64_t 
    //     int64_t seq = con.SequenceNumber();
    //     int64_t pos = cmesh->Block().Position(seq);
    //     int64_t posloc = pos - (nEquationsFull - nElements); // position in the local pressure solution
    //     mesh_sol(pos, 0) = pressure(posloc, 0);
    // }
    cmesh->TransferMultiphysicsSolution();
}

void PrintResults(TPZLinearAnalysis &an, TPZCompMesh *cmesh, int resolution)
{

    std::cout << "--------- Post Process ---------" << std::endl;
    TPZSimpleTimer postProc("Post processing time");
    const std::string plotfile = "postprocess";

    TPZVec<std::string> fields = {
        "Pressure",
        "ExactPressure",
        "Flux",
        "ExactFlux"};
    auto vtk = TPZVTKGenerator(cmesh, fields, plotfile, resolution);
    vtk.SetNThreads(global_nthread);
    vtk.Do();
    std::cout << "Total time = " << postProc.ReturnTimeDouble() / 1000. << " s" << std::endl;
}

void IterativeSolver(TPZGeoMesh* gmesh, ProblemData* inputData, REAL alpha, TPZAnalyticSolution* sol, std::ofstream &outfile)
{
    // velocity computational mesh
    TPZCompMesh *cmesh_v = CreateCMeshV(inputData, gmesh);

    // pressure computational mesh
    TPZCompMesh *cmesh_p = CreateCMeshP(inputData, gmesh);

    // multiphysics computational mesh
    TPZMultiphysicsCompMesh *cmesh_m = CreateMultiPhysicsMesh(inputData, gmesh, sol, true, alpha);
    
    // static condensation
    CondenseElements(inputData, cmesh_m, gmesh, true);
    {
        std::ofstream cmeshFile("cmesh_m.txt");
        cmesh_m->Print(cmeshFile);
    }

    // Number of equations without condense elements
    const int nEquationsFull = cmesh_m->Solution().Rows();
    std::cout << "Number of equations before condensation = = " << nEquationsFull << std::endl;
    if (!outfile.fail())
        outfile << "Number of equations before condensation = " << nEquationsFull << std::endl;

    //Number of condensed problem.
    int nEquationsCondensed = cmesh_m->NEquations();
    std::cout << "Number of equations after condensation = " << nEquationsCondensed << std::endl;
    if (!outfile.fail())
        outfile << "Number of equations after condensation = " << nEquationsCondensed << std::endl;

    //Create analysis environment
    TPZLinearAnalysis an(cmesh_m, RenumType::ENone);
    an.SetExact(sol->ExactSolution(),3);

    SolveProblemIterative(an, cmesh_m, alpha, 1.e-9, outfile);

    std::cout << "--------- PostProcess ---------" << std::endl;
    PrintResults(an, cmesh_m, inputData->Resolution());
}

void DirectSolver(TPZGeoMesh* gmesh, ProblemData* inputData, TPZAnalyticSolution* sol, std::ofstream &outfile)
{
    // velocity computational mesh
    TPZCompMesh *cmesh_v = CreateCMeshV(inputData, gmesh);

    // pressure computational mesh
    TPZCompMesh *cmesh_p = CreateCMeshP(inputData, gmesh);

    // multiphysics computational mesh
    TPZMultiphysicsCompMesh *cmesh_m = CreateMultiPhysicsMesh(inputData, gmesh, sol, false);
    
    // static condensation
    CondenseElements(inputData, cmesh_m, gmesh);

    // Number of equations without condense elements
    const int nEquationsFull = cmesh_m->Solution().Rows();
    std::cout << "Number of equations before condensation = = " << nEquationsFull << std::endl;
    if (!outfile.fail())
        outfile << "Number of equations before condensation = " << nEquationsFull << std::endl;

    //Number of condensed problem.
    int nEquationsCondensed = cmesh_m->NEquations();
    std::cout << "Number of equations after condensation = " << nEquationsCondensed << std::endl;
    if (!outfile.fail())
        outfile << "Number of equations after condensation = " << nEquationsCondensed << std::endl;

    //Create analysis environment
    TPZLinearAnalysis an(cmesh_m, RenumType::EMetis);
    an.SetExact(sol->ExactSolution(),3);

    SolveProblemDirect(an, cmesh_m, outfile);

    std::cout << "--------- PostProcess ---------" << std::endl;
    PrintResults(an, cmesh_m, inputData->Resolution());
}