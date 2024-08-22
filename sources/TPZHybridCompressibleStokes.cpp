#include <pzfmatrix.h>
#include <TPZBndCondT.h>
#include <pzaxestools.h>
#include <pzlog.h>
#include <fstream>
#include <ostream>
#include <iostream>
#include "TPZHybridCompressibleStokes.h"

#ifdef PZ_LOG
static TPZLogger logger("pz.stokesmaterial");
#endif

TPZHybridCompressibleStokes::TPZHybridCompressibleStokes() : TPZHybridStokes(), fAlpha(0.0)
{
}

TPZHybridCompressibleStokes::TPZHybridCompressibleStokes(int matID, int dimension, REAL viscosity, REAL alpha) : TPZHybridStokes(matID, dimension, viscosity), fAlpha(alpha)
{
}

TPZHybridCompressibleStokes::~TPZHybridCompressibleStokes()
{
}

void TPZHybridCompressibleStokes::Contribute(const TPZVec<TPZMaterialDataT<STATE>> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef)
{
    int64_t nShapeV = datavec[EVindex].fVecShapeIndex.NElements(); // number of velocity Hdiv shape functions
    TPZFNMatrix<150, REAL> PhiV(fDimension, nShapeV, 0.0);
    TPZFNMatrix<20, REAL> &divPhiV = datavec[EVindex].divphi;

    TPZFMatrix<REAL> &PhiP = datavec[EPindex].phi;
    int64_t nShapeP = PhiP.Rows(); // number of pressure H1 shape functions

    const int nterms = fDimension * (fDimension + 1) / 2;
    TPZFNMatrix<150, REAL> StrainRate(nterms, nShapeV, 0.0); // Using voight notation

    if (datavec[EVindex].fNeedsDeformedDirectionsFad)
    {
        for (int64_t j = 0; j < nShapeV; j++)
        {
            int cont = fDimension - 1;
            for (int64_t i = 0; i < fDimension; i++)
            {
                PhiV(i, j) = datavec[EVindex].fDeformedDirectionsFad(i, j).val();
                StrainRate(i, j) = datavec[EVindex].fDeformedDirectionsFad(i, j).fastAccessDx(i);
                for (int64_t k = i + 1; k < fDimension && k != i; k++)
                {
                    StrainRate(++cont, j) = 0.5 * (datavec[EVindex].fDeformedDirectionsFad(i, j).fastAccessDx(k) + datavec[EVindex].fDeformedDirectionsFad(k, j).fastAccessDx(i));
                }
            }
        }
    }

    TPZFNMatrix<9, REAL> VoightCorrection(nterms, nterms, 0.0); // This is to multiply by 2 the off diagonal part of strain rate tensor to account for its symmetry
    for (int64_t i = 0; i < nterms; i++)
    {
        VoightCorrection(i, i) = i < fDimension ? 1.0 : 2.0;
    }

    TPZFNMatrix<3, REAL> SourceTerm(fDimension, 1.0, 0.0);
    TPZVec<REAL> sourceAux(3);
    if (this->HasForcingFunction())
    {
        this->ForcingFunction()(datavec[EVindex].x, sourceAux);
        for (int64_t i = 0; i < fDimension; i++)
        {
            SourceTerm(i, 0) = sourceAux[i] * fViscosity;
        }
    }

    // Body Forces contribution
    ef.AddContribution(0, 0, PhiV, true, SourceTerm, false, weight);

    // Flux Matrix A contribution
    TPZFNMatrix<150, REAL> matrixC;
    VoightCorrection.Multiply(StrainRate, matrixC);

    REAL factor = 2.0 * fViscosity * weight;
    ek.AddContribution(0, 0, StrainRate, true, matrixC, false, factor);

    factor = fViscosity * weight;
    ek.AddContribution(0, 0, divPhiV, false, divPhiV, true, factor);

    // Divergence Matrix B contribution
    factor = -1.0 * weight;
    ek.AddContribution(0, nShapeV, divPhiV, false, PhiP, true, factor);

    // Divergence Matrix BT contribution
    ek.AddContribution(nShapeV, 0, PhiP, false, divPhiV, true, factor);

    // Compressibility matrix C contribution
    TPZFNMatrix<150, REAL> I(nShapeP,nShapeP,0.0);
    I.Identity();
    factor = -1.0 * fAlpha * weight;
    ek.AddContribution(nShapeV, nShapeV, I, false, I, false, factor);

    if (datavec.size() > 2)
    {
        TPZFMatrix<REAL> &phiG = datavec[EGindex].phi;
        TPZFMatrix<REAL> &phipM = datavec[EPMindex].phi;

        // Pressure and distributed flux
        for (int j = 0; j < nShapeP; j++)
        {
            ek(nShapeV + nShapeP, nShapeV + j) += PhiP(j, 0) * phiG(0, 0) * weight;
            ek(nShapeV + j, nShapeV + nShapeP) += PhiP(j, 0) * phiG(0, 0) * weight;
        }

        // Injection and average-pressure
        ek(nShapeV + nShapeP + 1, nShapeV + nShapeP) += phiG(0, 0) * phipM(0, 0) * weight;
        ek(nShapeV + nShapeP, nShapeV + nShapeP + 1) += phiG(0, 0) * phipM(0, 0) * weight;
    }

#ifdef PZ_LOG
    if (logger.isDebugEnabled())
    {
        std::stringstream sout;
        ek.Print("ek", sout, EMathematicaInput);
        ef.Print("ef", sout, EMathematicaInput);
        sout << std::endl
             << std::endl;
        LOGPZ_DEBUG(logger, sout.str())
    }
#endif
}