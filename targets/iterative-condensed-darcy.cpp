//
// Created by Giovane Avancini on 01/02/24.
//
#ifdef HAVE_CONFIG_H
#include <pz_config.h>
#endif

#include "meshpath_config.h"
#include <TPZGeoMeshTools.h>
#include "TPZAnalyticSolution.h"
#include <TPZGmshReader.h>
#include "TPZCompMeshTools.h"
#include "pzlog.h"
#include "pzshapequad.h"
#include "pzshapetriang.h"
#include "pzshapecube.h"
#include "pzshapetetra.h"
#include "TPZTimer.h"
#include "fstream"
#include "TPZSimpleTimer.h"
#include "TPZVTKGenerator.h"
#include "TPZRefPattern.h"
#include "tpzgeoelrefpattern.h"
#include "TPZRefPatternDataBase.h"
#include "pzbuildmultiphysicsmesh.h"
#include "TPZHDivApproxCreator.h"
#include "DarcyFlow/TPZMixedDarcyFlow.h"
#include <TPZSSpStructMatrix.h> //symmetric sparse matrix storage
#include <pzstepsolver.h> //for TPZStepSolver
#include "pzblockdiag.h"
#include "pzbdstrmatrix.h"
#include "pzcmesh.h"
#include "TPZMultiphysicsCompMesh.h"
#include "TPZLinearAnalysis.h"
#include <TPZVTKGeoMesh.h>
#include "ProblemData.h"
#include <TPZYSMPMatrix.h>
#include <TPZYSMPPardiso.h>
#include <TPZSYSMPPardiso.h>
#include "TPZMixedCompressibleDarcyFlow.h"
#include <pzcondensedcompel.h>
#include <pzelchdiv.h>

const int global_nthread = 8;
enum EMatid  {ENone, EDomain, EBoundary, EPont, EWrap, EIntface, EPressureHyb};

template<class tshape>
TPZGeoMesh* CreateGeoMesh(int xdiv, EMatid volId, EMatid bcId);

template<class tshape>
TPZGeoMesh* ReadMeshFromGmsh(std::string file_name);

void SolveProblemDirect(TPZLinearAnalysis &an, TPZCompMesh *cmesh);

void SolveProblemIterative(TPZLinearAnalysis &an, TPZCompMesh *cmesh, double alpha, double tol);

void PrintResults(TPZLinearAnalysis &an, TPZCompMesh *cmesh, int resolution = 0);

//Analytical solution
constexpr int solOrder{4};
auto exactSol = [](const TPZVec<REAL> &loc,
    TPZVec<STATE>&u,
    TPZFMatrix<STATE>&gradU){
    const auto &x=loc[0];
    const auto &y=loc[1];
    const auto &z=loc[2];

    REAL aux = 1./sinh(sqrt(2)*M_PI);
    u[0] = sin(M_PI*x)*sin(M_PI*y)*sinh(sqrt(2)*M_PI*z)*aux;
    gradU(0,0) = M_PI*cos(M_PI*x)*sin(M_PI*y)*sinh(sqrt(2)*M_PI*z)*aux;
    gradU(1,0) = M_PI*cos(M_PI*y)*sin(M_PI*x)*sinh(sqrt(2)*M_PI*z)*aux;
    gradU(2,0) = sqrt(2)*M_PI*cosh(sqrt(2)*M_PI*z)*sin(M_PI*x)*sin(M_PI*y)*aux;
};

int main(int argc, char *argv[])
{
    #ifdef PZ_LOG
    TPZLogger::InitializePZLOG(std::string(MESHES_DIR) + "/" + "log4cxx.cfg");
    #endif
    const int xdiv = 20;
    const int pOrder = 2;
    HDivFamily hdivfam = HDivFamily::EHDivConstant;
    
    // Creates/import a geometric mesh  
    TPZGeoMesh* gmesh = CreateGeoMesh<pzshape::TPZShapeCube>(xdiv, EDomain, EBoundary);
    int dim = gmesh->Dimension();
    TPZHDivApproxCreator hdivCreator(gmesh);
    hdivCreator.HdivFamily() = hdivfam;
    hdivCreator.ProbType() = ProblemType::EDarcy;
    hdivCreator.IsRigidBodySpaces() = false;
    hdivCreator.SetDefaultOrder(pOrder);
    hdivCreator.SetExtraInternalOrder(0);
    hdivCreator.SetShouldCondense(true);
    hdivCreator.SetShouldCondensePressure(true);
    hdivCreator.HybridType() = HybridizationType::ENone;
    // hdivCreator.HybridType() = HybridizationType::EStandard;

    // Prints gmesh mesh properties
    {
        std::string vtk_name = "geoMesh.vtk";
        std::ofstream vtkfile(vtk_name.c_str());
        TPZVTKGeoMesh::PrintGMeshVTK(gmesh, vtkfile, true);
    }

    //Insert Materials
    // TPZMixedDarcyFlow* matdarcy = new TPZMixedDarcyFlow(EDomain,dim);
    REAL alpha = 0.00001;
    TPZMixedCompressibleDarcyFlow* matdarcy = new TPZMixedCompressibleDarcyFlow(EDomain,dim,alpha);
    matdarcy->SetConstantPermeability(1.);
    matdarcy->SetExactSol(exactSol,4);

    hdivCreator.InsertMaterialObject(matdarcy);

    TPZFMatrix<STATE> val1(1,1,0.);
    TPZManVector<STATE> val2(1,0.);
    TPZBndCondT<STATE> *BCond1 = matdarcy->CreateBC(matdarcy, EBoundary, 0, val1, val2);
    BCond1->SetForcingFunctionBC(exactSol,4);
    hdivCreator.InsertMaterialObject(BCond1);

    //Multiphysics mesh
    TPZMultiphysicsCompMesh *cmesh = hdivCreator.CreateApproximationSpace();
    {
        std::string txt = "cmesh.txt";
        std::ofstream myfile(txt);
        cmesh->Print(myfile);
    }
  
    // Number of equations without condense elements
    const int nEquationsFull = cmesh->Solution().Rows();
    std::cout << "Number of equations = " << nEquationsFull << std::endl;

    TPZTimer clock,clock2;
    clock.start();

    //Number of condensed problem.
    int nEquationsCondensed = cmesh->NEquations();
    std::cout << "Number of equations condensed = " << nEquationsCondensed << std::endl;
    //Create analysis environment
    TPZLinearAnalysis an(cmesh, RenumType::EMetis);
    an.SetExact(exactSol,solOrder);

    bool domHyb = false;
    // SolveProblemDirect(an, cmesh);
    SolveProblemIterative(an, cmesh, alpha, 1.e-9);
    clock.stop();

    std::cout << "--------- PostProcess ---------" << std::endl;
    PrintResults(an, cmesh, 0);

    return 0;
}

template <class tshape>
TPZGeoMesh* CreateGeoMesh(int xdiv, EMatid volId, EMatid bcId)
{
    MMeshType meshType;
    int dim = tshape::Dimension;
    TPZVec<int> nDivs;

    if (dim == 2) nDivs = {xdiv,xdiv};
    if (dim == 3) nDivs = {xdiv,xdiv,xdiv};

    switch (tshape::Type())
    {
    case ETriangle:
        meshType = MMeshType::ETriangular;
        break;
    case EQuadrilateral:
        meshType = MMeshType::EQuadrilateral;
        break;
    case ETetraedro:
        meshType = MMeshType::ETetrahedral;
        break;
    case ECube:
        meshType = MMeshType::EHexahedral;
        break;
        case EPrisma:
        meshType = MMeshType::EPrismatic;
        break;
    default:
        DebugStop();
    }

    TPZManVector<REAL,3> minX = {0,0,0};
    TPZManVector<REAL,3> maxX = {1,1,1};
    int nMats = 2*dim+1;

    //all bcs share the same id
    constexpr bool createBoundEls{true};
    TPZVec<int> matIds(nMats,bcId);
    matIds[0] = volId;
    
    TPZGeoMesh* gmesh = TPZGeoMeshTools::CreateGeoMeshOnGrid(dim, minX, maxX, matIds, nDivs, meshType,createBoundEls);
    
    return gmesh;
}

void SolveProblemDirect(TPZLinearAnalysis &an, TPZCompMesh *cmesh)
{
    TPZSSpStructMatrix<STATE,TPZStructMatrixOR<STATE>> matskl(cmesh);   
    
    matskl.SetNumThreads(global_nthread);
    
    an.SetStructuralMatrix(matskl);
    
    ///Setting a direct solver
    TPZStepSolver<STATE> step;
    step.SetDirect(ELDLt);//ELU //ECholesky // ELDLt
    
    an.SetSolver(step);

    //assembles the system
    std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
    an.Assemble();
    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
    std::cout << "Time Assemble = " << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << "[ms]" << std::endl;

    ///solves the system
    std::chrono::steady_clock::time_point begin2 = std::chrono::steady_clock::now();
    an.Solve();
    std::chrono::steady_clock::time_point end2 = std::chrono::steady_clock::now();
    std::cout << "Time Solve = " << std::chrono::duration_cast<std::chrono::milliseconds>(end2 - begin2).count() << "[ms]" << std::endl;
}

void SolveProblemIterative(TPZLinearAnalysis &an, TPZCompMesh *cmesh, double alpha, double tol)
{
    TPZSSpStructMatrix<STATE,TPZStructMatrixOR<STATE>> matskl(cmesh);   
    
    matskl.SetNumThreads(global_nthread);
    
    an.SetStructuralMatrix(matskl);
    
    ///Setting a direct solver
    TPZStepSolver<STATE> step;
    step.SetDirect(ECholesky);//ELU //ECholesky // ELDLt
    
    an.SetSolver(step);

    //assembles the system
    std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
    an.Assemble();
    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
    std::cout << "Time Assemble = " << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << "[ms]" << std::endl;

    //Computing the number of hdiv domain elements, the number of facet connects and the number of facet equations
    TPZMultiphysicsCompMesh* cmesh_m = dynamic_cast<TPZMultiphysicsCompMesh*>(cmesh);
    TPZCompMesh* cmesh_u = cmesh_m->MeshVector()[0];
    int64_t nel_u = cmesh_u->NElements();
    int64_t nElements = 0;
    int64_t nFacetConnects = 0;
    int64_t nFacetEqs = 0;

    for (int64_t iel = 0; iel < nel_u; iel++)
    {
        TPZCompEl *cel = cmesh_u->Element(iel);
        if (!cel) continue;
        if (cel->Dimension() != cmesh->Dimension()) continue;
        
        TPZGeoEl *gel = cel->Reference();
        nFacetConnects += gel->NSides(cmesh->Dimension()-1);
        nElements++;
    }
    
    for (auto& c : cmesh_u->ConnectVec())
    {
        if (c.LagrangeMultiplier() != 0 || c.HasDependency() || c.IsCondensed() || c.NElConnected() == 1) continue;
        nFacetEqs += c.NShape();
    }

    /*
    For Hdiv constant space, the global matrix B has dimensions (nFacetEqs x nElements)
    and only has contribution of the constant flux dof of each face equal to -fsideorient.
    Thus the number of nonzeros of B is nFaceConnects.
    */
    TPZVec<int64_t> iBT(nElements+1,0);
    TPZVec<int64_t> jBT(nFacetConnects,0);
    TPZVec<STATE> valBT(nFacetConnects,0.);

    int64_t row = 0;
    int count = 0;
    for (int64_t iel = 0; iel < nel_u; iel++)
    {
        int64_t col = 0;
        TPZCompEl *cel = cmesh_u->Element(iel);
        if (!cel) continue;
        if (cel->Dimension() != cmesh->Dimension()) continue;
        TPZGeoEl *gel = cel->Reference();
        int nsides = gel->NSides();
        int nfacets = gel->NSides(gel->Dimension()-1);
        int eqcont = 0;
        iBT[row+1] = iBT[row] + nfacets;
        std::set<int64_t> connectlist;
        cel->BuildConnectList(connectlist);
        for (auto it = connectlist.begin(); it != connectlist.end(); it++)
        {
            TPZConnect &con = cmesh_m->ConnectVec()[*it];

            if (con.LagrangeMultiplier() != 0 || con.HasDependency() || con.IsCondensed() || con.NElConnected() == 1) continue;

            for (int iside = gel->NCornerNodes(); iside < nsides-1; iside++)
            {
                if (gel->SideDimension(iside) != cmesh->Dimension() - 1) continue; // only the facets
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
        row++;
    }

    TPZFYsmpMatrixPardiso<STATE> BT(nElements,nFacetEqs);
    BT.SetData(iBT,jBT,valBT);
    
    //Gettig the Global stiffness matrix from solver
    auto KG = an.MatrixSolver<STATE>().Matrix();
    
    TPZFMatrix<STATE> force = an.Rhs();
    TPZFMatrix<STATE> rhs(cmesh->NEquations(),1,0.);
    TPZFMatrix<STATE> dsol(cmesh->NEquations(),1,0.);

    //Obtaining the initial solution
    begin = std::chrono::steady_clock::now();
    an.Solve();
    TPZFMatrix<STATE> sol = an.Solution();
    end = std::chrono::steady_clock::now();
    std::cout << "Time for Initial Solver = " << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << "[ms]" << std::endl;
    
    //Computing the residual (force is not used)
    TPZFMatrix<STATE> aux(nElements,1,0.);
    BT.MultAdd(sol, force, aux, 1.0, 0.0); //aux = BT * sol
    BT.MultAdd(aux, force, rhs, -1.0/alpha, 0.0, 1); //rhs = -Transpose(BT) * aux / alpha = -Transpose(BT) * BT * sol / alpha
    dsol = rhs;

    //Computing initial pressure (force is not used)
    TPZFMatrix<STATE> pressure(nElements,1,0.);
    BT.MultAdd(sol, force, pressure, 1.0/alpha, 0.0); //pressure = BT * sol / alpha
    
    double norm_dsol=1.0, norm_rhs=1.0;
    int nit = 0;
    const int size = rhs.Rows();
    begin = std::chrono::steady_clock::now();
    while (norm_dsol > tol || norm_rhs > tol)
    {
        KG->SolveDirect(dsol,ECholesky);
        sol += dsol;
        //Updating the residual (force is not used)
        BT.MultAdd(sol, force, aux, 1.0, 0.0); //aux = BT * sol
        BT.MultAdd(aux, force, rhs, -1.0/alpha, 0.0, 1); //rhs = Transpose(BT) * aux = Transpose(BT) * BT * sol
        //Updating pressure (force is not used)
        TPZFMatrix<STATE> dp(nElements,1,0.);
        BT.MultAdd(sol, force, dp, 1.0/alpha, 0.0); //dp = BT * sol / alpha
        pressure += dp;
        nit++;
        norm_dsol = 0.0;
        norm_rhs = 0.0;
        for (int64_t i = 0; i < size; i++)
        {
            norm_dsol += dsol(i,0) * dsol(i,0);
            norm_rhs += rhs(i,0) * rhs(i,0);
        }
        norm_dsol = sqrt(norm_dsol);
        norm_rhs = sqrt(norm_rhs);
        dsol=rhs;
        std::cout << "Iteration: " << nit << ". dsol_norm: " << norm_dsol << ", rhs_norm: " << norm_rhs << std::endl;

        if (nit > 200)
        {
            std::cout << "Solver diverged.\n";
            break;
        }
    }
    
    end = std::chrono::steady_clock::now();
    std::cout << "Time Iterative Process = " << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << "[ms]" << std::endl;

    //Transfering the solution to the mesh
    an.Solution() = sol;
    an.LoadSolution();

    TPZFMatrix<REAL>& mesh_sol = cmesh->Solution();
    const int64_t nEquationsFull = mesh_sol.Rows();
    for (const TPZConnect& con : cmesh->ConnectVec())
    {
        if (con.LagrangeMultiplier() != 1) continue; //only pressure connect
        int64_t seq = con.SequenceNumber();
        int64_t pos = cmesh->Block().Position(seq);
        int64_t posloc = pos - (nEquationsFull - nElements); //position in the local pressure solution
        mesh_sol(pos,0) = pressure(posloc, 0);
    }
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