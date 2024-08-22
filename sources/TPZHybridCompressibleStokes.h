#include <pzfmatrix.h>
#include "TPZHybridStokes.h"
#include <math.h>

#ifndef TPZHYBRIDCOMPRESSIBLESTOKES
#define TPZHYBRIDCOMPRESSIBLESTOKES

class TPZHybridCompressibleStokes : public TPZHybridStokes{
    
protected:
    /// compressibility parameter
    REAL fAlpha;
    
public:
    /// Empty Constructor
    TPZHybridCompressibleStokes();

    /// Creates a material object and inserts it in the vector of material pointers of the mesh
    TPZHybridCompressibleStokes(int matID, int dimension, REAL viscosity, REAL alpha);
    
    /// Destructor
    ~TPZHybridCompressibleStokes();
    
    // Contribute Methods being used - Multiphysics
    virtual void Contribute(const TPZVec<TPZMaterialDataT<STATE>> &datavec,
                            REAL weight,TPZFMatrix<STATE> &ek,
                            TPZFMatrix<STATE> &ef) override;
};

#endif