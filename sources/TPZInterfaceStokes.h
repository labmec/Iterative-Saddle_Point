#include <pzfmatrix.h>
#include <TPZMaterial.h>
#include <TPZMatBase.h>
#include <TPZMaterialData.h>
#include <TPZMaterialDataT.h>
#include <TPZMatCombinedSpaces.h>
#include <TPZMatInterfaceCombinedSpaces.h>
#include <TPZLagrangeMultiplier.h>
#include <math.h>

#ifndef TPZINTERFACEMATERIAL_H
#define TPZINTERFACEMATERIAL_H

class TPZInterfaceStokes : public TPZMatBase<STATE, TPZMatCombinedSpacesT<STATE>, TPZMatInterfaceCombinedSpaces<STATE>>,
    
    public TPZLagrangeMultiplierBase
    
{
    
    using TBase = TPZMatBase<STATE, TPZMatCombinedSpacesT<STATE>,TPZMatInterfaceCombinedSpaces<STATE>>;
    
    bool IsLagrangeMult() override{return true;};
    
protected:
    int fDimension;
    int fVindex = 0;
    int fPindex = 1;
    int fNStateVariables = -1;
    STATE fMultiplier = -1;
    
    REAL fBigNumber = pow(10,std::numeric_limits<STATE>::max_digits10*2/3);
    bool fIsAxisymmetric;

public:
    /// Creates a material object
    TPZInterfaceStokes(int matID, int dimension, bool isAxisymmetric = false);
    
    /// Destructor
    ~TPZInterfaceStokes();
    
    STATE InnerProductVec(TPZFMatrix<STATE>& S, TPZFMatrix<STATE>& T);
    
    int Dimension() const override {return fDimension;}
    
    void ContributeInterface(const TPZMaterialDataT<STATE>& data,
                             const std::map<int, TPZMaterialDataT<STATE>>& dataleft,
                             const std::map<int, TPZMaterialDataT<STATE>>& dataright, REAL
                             weight,
                             TPZFMatrix<STATE>& ek,
                             TPZFMatrix<STATE>& ef) override;
    
    void ContributeBCInterface(const TPZMaterialDataT<STATE> &data,
                               const std::map<int, TPZMaterialDataT<STATE>> &dataleft,
                               REAL weight,
                               TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef,
                               TPZBndCondT<STATE> &bc) override;
    
    void FillDataRequirementsInterface(TPZMaterialDataT<STATE> &data,
                                       std::map<int, TPZMaterialDataT<STATE>> &datavec_left,
                                       std::map<int, TPZMaterialDataT<STATE>> &datavec_right) override;
    
    void Contribute(const TPZVec<TPZMaterialDataT<STATE>> &datavec,
                    REAL weight,
                    TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef) override;

    void ContributeBC(const TPZVec<TPZMaterialDataT<STATE>> &datavec,
                      REAL weight, TPZFMatrix<STATE> &ek,
                      TPZFMatrix<STATE> &ef,
                      TPZBndCondT<STATE> &bc) override;
    
    void SolutionInterface(const TPZMaterialDataT<STATE> &data,
                           const std::map<int, TPZMaterialDataT<STATE>> &dataleftvec,
                           const std::map<int, TPZMaterialDataT<STATE>> &datarightvec,
                           int var, TPZVec<STATE> &Solout) override;


    void SolutionInterface(const TPZMaterialDataT<STATE> &data,
                           const std::map<int, TPZMaterialDataT<STATE>> &dataleftvec,
                           const std::map<int, TPZMaterialDataT<STATE>> &datarightvec,
                           int var, TPZVec<STATE> &Solout,
                           TPZCompEl *left,TPZCompEl *right) override;
    
    void Solution(const TPZVec<TPZMaterialDataT<STATE>> &datavec,
                  int var, TPZVec<STATE> &sol) override;

    int GetIntegrationOrder(const TPZVec<int> &porder_left, const TPZVec<int> &porder_right) const override;

    int NStateVariables() const override
    {return fNStateVariables;}
    
    virtual void SetMultiplier(STATE mult){
        fMultiplier = mult;
    }
    
};
#endif