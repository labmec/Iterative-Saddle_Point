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
#include <TPZSYSMPPardiso.h>
#include "TPZMixedCompressibleDarcyFlow.h"
#include <pzcondensedcompel.h>
#include <pzelchdiv.h>

const int global_nthread = 0;
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
    const int xdiv = 2;
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
    // hdivCreator.SetShouldCondense(false);
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
    TPZMixedCompressibleDarcyFlow* matdarcy = new TPZMixedCompressibleDarcyFlow(EDomain,dim,0.0);
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
    const int nEquationsFull = cmesh->NEquations();
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
    SolveProblemIterative(an, cmesh, 0.01, 1.e-9);
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

    auto K = an.MatrixSolver<STATE>().Matrix();
    {
        std::ofstream file("matrix");
        K->Print("GK", file, EMathematicaInput);
    }

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
    step.SetDirect(ELDLt);//ELU //ECholesky // ELDLt
    
    an.SetSolver(step);

    //assembles the system
    std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
    an.Assemble();
    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
    std::cout << "Time Assemble = " << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << "[ms]" << std::endl;

    //Computing the number of hdiv domain elements, the number of face connects and their number of hdiv face functions
    TPZMultiphysicsCompMesh* cmesh_m = dynamic_cast<TPZMultiphysicsCompMesh*>(cmesh);
    TPZCompMesh* cmesh_u = cmesh_m->MeshVector()[0];
    int64_t nel_u = cmesh_u->NElements();
    int64_t nElements = 0;
    int64_t nFaceConnects = 0;
    int64_t maxElementNeq = 0;

    for (int64_t iel = 0; iel < nel_u; iel++)
    {
        TPZCompEl *cel = cmesh_u->Element(iel);
        if (!cel) continue;
        if (cel->Dimension() != cmesh->Dimension()) continue;
        
        TPZGeoEl *gel = cel->Reference();
        nFaceConnects += gel->NSides(cmesh->Dimension()-1);
        nElements++;
        int nconnects = cel->NConnects();
        int nElementShape = 0;
        for (int ic = 0; ic < nconnects; ic++)
        {
            TPZConnect con = cel->Connect(ic);
            if (con.IsCondensed() || con.HasDependency() || con.NElConnected() == 1) continue;
            nElementShape += cel->Connect(ic).NShape();
        }
        if (nElementShape > maxElementNeq)
        {
            maxElementNeq = nElementShape;
        }
    }

    /*
    For Hdiv constant space, the global matrix B has dimensions (maxElementNeq x nElements)
    and only has contribution of the constant flux dof of each face equal to -fsideorient.
    Thus the number of nonzeros of B is nFaceConnects.
    */
    TPZVec<int64_t> iB(nElements+1,0);
    TPZVec<int64_t> jB(nFaceConnects,0);
    TPZVec<STATE> valB(nFaceConnects,0.);

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
        iB[row+1] = iB[row] + nfacets;
        
        for (int iside = gel->NCornerNodes(); iside < nsides-1; iside++)
        {
            if (gel->SideDimension(iside) != cmesh->Dimension() - 1) continue; //only the facets
            REAL side_orient = gel->NormalOrientation(iside);
            valB[count] = -side_orient;
            jB[count] = col;

            //Sum to the row index the number of shape functions of the connect
            TPZCompElSide celside(cel,iside);
            int connect_id = celside.ConnectIndex();
            TPZConnect con = cmesh_u->ConnectVec()[connect_id];
            col += con.NShape();
            count++; 
        }
        row++;
    }
    TPZFYsmpMatrix<STATE> B(nElements,maxElementNeq);
    B.SetData(iB,jB,valB);
    {
        std::ofstream file("matrix_B.txt");
        B.Print("B", file, EMathematicaInput);
    }
    
    //Gettig the Global stiffness matrix from solver
    auto KG = an.MatrixSolver<STATE>().Matrix();
    {
        std::ofstream file("matrix_original.txt");
        KG->Print("GK", file, EMathematicaInput);
        an.Rhs().Print("rhs", file, EMathematicaInput);
    }

    //Copying KG to K_it
    auto* K_it = dynamic_cast<TPZSYsmpMatrixPardiso<STATE>*>(KG->Clone());
    if (!K_it) DebugStop();

    //Adding -alpha to the zeroes diagonals of K_it
    for (TPZConnect con : cmesh->ConnectVec())
    {
        if (con.LagrangeMultiplier() != 1 || con.IsCondensed()) continue;
        int64_t seq = con.SequenceNumber();
        int64_t eq = cmesh->Block().Position(seq);
        K_it->PutVal(eq, eq, -alpha);
    }
    {
        std::ofstream file("matrix_compressible.txt");
        K_it->Print("Kcomp", file, EMathematicaInput);
    }
    
    TPZFMatrix<STATE> force = an.Rhs();
    TPZFMatrix<STATE> rhs = force;
    TPZFMatrix<STATE> dsol = rhs;
    TPZFMatrix<STATE> sol(rhs.Rows(),1,0.);

    //Solving
    double norm_dsol=1.0, norm_rhs=1.0;
    int nit = 0;
    const int size = rhs.Rows();
    while (norm_dsol > tol || norm_rhs > tol)
    {
        K_it->SolveDirect(dsol,ELDLt);
        sol += dsol;
        KG->MultAdd(sol, force, rhs, -1.0, 1.0); //Computing the residual rhs = F - KG.sol

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

        if (nit > 50)
        {
            std::cout << "Solver diverged.\n";
            break;
        }
    }
    
    //Transfering the solution to the mesh
    an.Solution() = sol;
    an.LoadSolution();
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