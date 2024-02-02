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

const int global_nthread = 0;
enum EMatid  {ENone, EDomain, EBoundary, EPont, EWrap, EIntface, EPressureHyb};

template<class tshape>
TPZGeoMesh* CreateGeoMesh(int xdiv, EMatid volId, EMatid bcId);

template<class tshape>
TPZGeoMesh* ReadMeshFromGmsh(std::string file_name);

void SolveProblemDirect(TPZLinearAnalysis &an, TPZCompMesh *cmesh, bool filterEquations, bool domainHybridization);

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
    const int xdiv = 1;
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
    hdivCreator.HybridType() = HybridizationType::ESemi;
    // hdivCreator.HybridType() = HybridizationType::EStandard;

    // Prints gmesh mesh properties
    {
        std::string vtk_name = "geoMesh.vtk";
        std::ofstream vtkfile(vtk_name.c_str());
        TPZVTKGeoMesh::PrintGMeshVTK(gmesh, vtkfile, true);
    }

    //Insert Materials
    TPZMixedDarcyFlow* matdarcy = new TPZMixedDarcyFlow(EDomain,dim);
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
    SolveProblemDirect(an, cmesh, false, domHyb);
    clock.stop();

    std::cout << "--------- PostProcess ---------" << std::endl;
    PrintResults(an, cmesh, 2);

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

void SolveProblemDirect(TPZLinearAnalysis &an, TPZCompMesh *cmesh, bool filterEquations, bool domainHybridization)
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