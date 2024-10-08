#include <pzfmatrix.h>
#include <TPZBndCondT.h>
#include <pzaxestools.h>
#include <pzlog.h>
#include <fstream>
#include <ostream>
#include <iostream>

#include "TPZHybridStokes.h"

#ifdef PZ_LOG
static TPZLogger logger("pz.stokesmaterial");
#endif

TPZHybridStokes::TPZHybridStokes() : TBase()
{
}

TPZHybridStokes::TPZHybridStokes(int matID, int dimension, REAL viscosity) : TBase(matID), fDimension(dimension), fViscosity(viscosity)
{
}

TPZHybridStokes::~TPZHybridStokes()
{
}

void TPZHybridStokes::Contribute(const TPZVec<TPZMaterialDataT<STATE>> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef)
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

void TPZHybridStokes::ContributeBC(const TPZVec<TPZMaterialDataT<STATE>> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCondT<STATE> &bc)
{

    TPZFNMatrix<150, REAL> PhiV = datavec[EVindex].phi;
    TPZFMatrix<REAL> &PhiP = datavec[EPindex].phi;

    int64_t nShapeV = PhiV.Rows();
    int64_t nShapeP = PhiP.Rows();

    TPZFNMatrix<20, STATE> val1(3, 3, 0.); // grad
    TPZManVector<STATE, 4> val2(3, 0.);    // value
    REAL pressure = 0.0;

    if (bc.HasForcingFunctionBC())
    {
        TPZVec<STATE> vVal;
        TPZFMatrix<STATE> gradval;
        bc.ForcingFunctionBC()(datavec[EVindex].x, val2, val1);
        pressure = val2[3];
    }
    else
    {
        val1 = bc.Val1();
        val2 = bc.Val2();
    }

    switch (bc.Type())
    {
    case ENormVelocity: // Normal Velocity
    {
        REAL v_n = val2[0];            // if bc was set in .json file, the normal value was already pescribed
        if (bc.HasForcingFunctionBC()) // if the bc is set through an analytic solution, we need to compute its normal component
        {
            v_n = 0.0;
            for (int i = 0; i < fDimension; i++)
            {
                v_n += val2[i] * datavec[EVindex].normal[i];
            }
        }

        REAL factor = fBigNumber * weight;

        for (int64_t j = 0; j < nShapeV; j++)
        {
            ef(j) += v_n * PhiV(j, 0) * factor;

            for (int64_t i = 0; i < nShapeV; i++)
            {
                ek(i, j) += PhiV(i, 0) * PhiV(j, 0) * factor;
            }
        }
        break;
    }

    case ETanVelocity: // Tangential Velocity
    {
        TPZManVector<REAL, 3> v_t = {0., 0., 0.}; // for tangential bc, a vector is prescribed, so we take the inner product with the local tangential axe
        for (int i = 0; i < fDimension - 1; i++)
        {
            for (int j = 0; j < fDimension; j++)
            {
                v_t[i] += val2[j] * datavec[EVindex].axes(i, j);
            }
        }

        REAL factor = fBigNumber * weight;

        for (int64_t j = 0; j < nShapeV; j++)
        {
            for (int64_t k = 0; k < fDimension - 1; k++)
            {
                int64_t index1 = (fDimension - 1) * j + k;
                ef(index1) += -v_t[k] * PhiV(j, 0) * factor;

                for (int64_t i = 0; i < nShapeV; i++)
                {
                    for (int64_t l = 0; l < fDimension - 1; l++)
                    {
                        int64_t index2 = (fDimension - 1) * i + l;
                        if (k != l)
                            continue;
                        ek(index1, index2) += PhiV(i, 0) * PhiV(j, 0) * factor;
                    }
                }
            }
        }
        break;
    }

    case ENormStress: // Normal Stress
    {
        REAL sigma_nn = val2[0]; // if bc was set in .json file, the normal value was already prescribed
        if (bc.HasForcingFunctionBC())
        {                                                    // if the bc is set through an analytic solution, we need to compute its normal component from the velocity gradient
            const int n = fDimension * (fDimension + 1) / 2; // size of the vector in Voight notation

            TPZFNMatrix<6, REAL> sigmaVoight(n, 1, 0.);
            TPZFNMatrix<9, STATE> sigma(3, 3, 0.0);

            StressTensor(val1, sigmaVoight, pressure);
            FromVoigt(sigmaVoight, sigma);

            TPZFNMatrix<3, REAL> sigma_n(fDimension, 1, 0.0);

            for (int i = 0; i < fDimension; i++)
            {
                for (int j = 0; j < fDimension; j++)
                {
                    sigma_n(i, 0) += sigma(i, j) * datavec[EVindex].normal[j];
                }
            }

            sigma_nn = 0.0;
            for (int i = 0; i < fDimension; i++)
            {
                sigma_nn += sigma_n[i] * datavec[EVindex].normal[i];
            }
        }

        for (int64_t i = 0; i < nShapeV; i++)
        {
            ef(i) += sigma_nn * PhiV(i, 0) * weight;
        }
        break;
    }

    case ETanStress: // Tangential Stress
    {
        TPZManVector<REAL, 3> sigma_nt(fDimension - 1, 0.);
        if (bc.HasForcingFunctionBC())
        {
            const int n = fDimension * (fDimension - 1) / 2;

            TPZFNMatrix<6, REAL> sigmaVoight(n, 1, 0.0);
            TPZFNMatrix<9, STATE> sigma(3, 3, 0.0);

            StressTensor(val1, sigmaVoight, pressure);
            FromVoigt(sigmaVoight, sigma);

            TPZManVector<REAL, 3> sigma_n(fDimension, 0.0);

            for (int i = 0; i < fDimension; i++)
            {
                for (int j = 0; j < fDimension; j++)
                {
                    sigma_n[i] += sigma(i, j) * datavec[EVindex].normal[j];
                }
            }

            for (int i = 0; i < fDimension - 1; i++)
            {
                for (int j = 0; j < fDimension; j++)
                {
                    sigma_nt[i] += sigma_n[i] * datavec[EVindex].axes(i, j);
                }
            }
        }
        else
        {
            for (int i = 0; i < fDimension - 1; i++)
                for (int j = 0; j < fDimension; j++)
                    sigma_nt[i] += val2[j] * datavec[EVindex].axes(i, j);
        }

        for (int64_t i = 0; i < nShapeV; i++)
        {
            for (int j = 0; j < fDimension - 1; j++)
            {
                int64_t index = (fDimension - 1) * i + j;
                ef(index) += -PhiV(i, 0) * sigma_nt[j] * weight;
            }
        }

        break;
    }

    default:
    {
        std::cout << "ERROR: BOUNDARY NOT IMPLEMENTED" << std::endl;
        DebugStop();
        break;
    }
    }
}

int TPZHybridStokes::VariableIndex(const std::string &name) const
{

    if (!strcmp("Pressure", name.c_str()))
        return EPressure;
    if (!strcmp("ExactPressure", name.c_str()))
        return EExactPressure;
    if (!strcmp("ErrorPressure", name.c_str()))
        return EErrorPressure;

    if (!strcmp("Velocity", name.c_str()))
        return EVelocity;
    if (!strcmp("Flux", name.c_str()))
        return EVelocity;
    if (!strcmp("ExactVelocity", name.c_str()))
        return EExactVelocity;
    if (!strcmp("ExactFlux", name.c_str()))
        return EExactVelocity;
    if (!strcmp("ErrorVelocity", name.c_str()))
        return EErrorVelocity;
    if (!strcmp("ErrorFlux", name.c_str()))
        return EErrorVelocity;

    if (!strcmp("SourceTerm", name.c_str()))
        return ESourceTerm;

    if (!strcmp("Stress", name.c_str()))
        return EStress;
    if (!strcmp("ExactStress", name.c_str()))
        return EExactStress;
    if (!strcmp("ErrorStress", name.c_str()))
        return EErrorStress;

    std::cout << "\n\nVar index not implemented\n\n";
    DebugStop();

    return 0;
}

int TPZHybridStokes::NSolutionVariables(int var) const
{

    int aux;

    switch (var)
    {
    case EPressure:
    case EExactPressure:
    case EErrorPressure:
        aux = 1;
        break;

    case EVelocity:
    case EExactVelocity:
    case EErrorVelocity:
    case ESourceTerm:
        aux = 3;
        break;

    case EStress:
    case EExactStress:
    case EErrorStress:
        aux = 9;
        break;
    default:
        std::cout << "\n\nVar index not implemented!!!\n\n";
        DebugStop();
    }
    return aux;
}

void TPZHybridStokes::Solution(const TPZVec<TPZMaterialDataT<STATE>> &datavec, int var, TPZVec<STATE> &Solout)
{

    const int n = fDimension * (fDimension + 1) / 2;

    TPZManVector<STATE, 3> u_h = datavec[EVindex].sol[0];
    TPZFNMatrix<10, STATE> gradU_h = datavec[EVindex].dsol[0];

    TPZManVector<STATE, 3> p_h = datavec[EPindex].sol[0];

    TPZManVector<STATE, 4> sol_exact(4);
    TPZFNMatrix<9, STATE> gradsol_exact(3, 3);
    STATE p_exact = 0.;

    if (this->HasExactSol())
    {
        fExactSol(datavec[EVindex].x, sol_exact, gradsol_exact);
        p_exact = sol_exact[3];
    }

    Solout.Resize(NSolutionVariables(var));

    switch (var)
    {

    case EPressure: // pressure
    {
        Solout[0] = p_h[0];
    }
    break;

    case EExactPressure: // exact pressure
    {
        Solout[0] = p_exact;
    }
    break;

    case EErrorPressure: // pressure element error
    {
        Solout[0] = abs(p_exact - p_h[0]);
    }
    break;

    case EVelocity: // velocity
    {
        Solout[0] = u_h[0]; // vx
        Solout[1] = u_h[1]; // vy
        Solout[2] = u_h[2]; // vz
    }
    break;

    case EExactVelocity: // exact velocity
    {
        Solout[0] = sol_exact[0];
        Solout[1] = sol_exact[1];
        Solout[2] = sol_exact[2];
    }
    break;

    case EErrorVelocity: // velocity element error
    {
        Solout[0] = abs(sol_exact[0] - u_h[0]); // vx
        Solout[1] = abs(sol_exact[1] - u_h[1]); // vy
        Solout[2] = abs(sol_exact[2] - u_h[2]); // vz
    }
    break;

    case ESourceTerm: // source term
    {
        TPZManVector<STATE, 3> f(3, 0.);

        if (!this->HasForcingFunction())
            this->ForcingFunction()(datavec[EVindex].x, f);

        Solout[0] = f[0]; // fx
        Solout[1] = f[1]; // fy
        Solout[2] = f[2]; // fz
    }
    break;

    case EStress: // stress
    {
        TPZFNMatrix<6, STATE> sigmaVoigt(n, 1, 0.);
        TPZFNMatrix<9, STATE> sigma(3, 3, 0.0);

        StressTensor(gradU_h, sigmaVoigt, p_h[0]);
        FromVoigt(sigmaVoigt, sigma);

        for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++)
                Solout[i * 3 + j] = sigma(i, j);
    }
    break;

    case EExactStress: // exact stress
    {
        TPZFNMatrix<6, STATE> sigmaVoigt_exact(n, 1, 0.0);
        TPZFNMatrix<9, STATE> sigma_exact(3, 3, 0.0);

        StressTensor(gradsol_exact, sigmaVoigt_exact, p_exact);
        FromVoigt(sigmaVoigt_exact, sigma_exact);

        for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++)
                Solout[i * 3 + j] = sigma_exact(i, j);
    }
    break;

    case EErrorStress: // stress element error
    {
        TPZFNMatrix<6, STATE> sigmaVoigt(n, 1, 0.0), sigmaVoigt_exact(n, 1, 0.0);
        TPZFNMatrix<9, STATE> sigma_h(3, 3, 0.0), sigma_exact(3, 3, 0.0);

        StressTensor(gradU_h, sigmaVoigt, p_h[0]);
        StressTensor(gradsol_exact, sigmaVoigt_exact, p_exact);

        FromVoigt(sigmaVoigt, sigma_h);
        FromVoigt(sigmaVoigt_exact, sigma_exact);

        for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++)
                Solout[i * 3 + j] = abs(sigma_exact(i, j) - sigma_h(i, j));
    }
    break;

    default:
    {
        std::cout << "\n\nVar index not implemented\n\n";
        DebugStop();
    }
    }
}

void TPZHybridStokes::FillDataRequirements(TPZVec<TPZMaterialDataT<STATE>> &datavec) const
{
    int64_t ndata = datavec.size();
    for (int idata = 0; idata < ndata; idata++)
    {
        datavec[idata].SetAllRequirements(false);
        datavec[idata].fNeedsSol = true;
        datavec[idata].fNeedsHSize = true;
        datavec[idata].fNeedsNormal = true;
    }
    datavec[0].fNeedsDeformedDirectionsFad = true;
}

void TPZHybridStokes::FillBoundaryConditionDataRequirements(int type, TPZVec<TPZMaterialDataT<STATE>> &datavec) const
{
    datavec[EVindex].fNeedsSol = false;
    datavec[EPindex].fNeedsSol = false;
    datavec[EVindex].fNeedsNormal = true;
    datavec[EPindex].fNeedsNormal = true;
}

void TPZHybridStokes::Errors(const TPZVec<TPZMaterialDataT<STATE>> &data, TPZVec<REAL> &errors)
{
    // 0: L2 p, 1: L2 p_ex, 2: L2 u, 3: L2 u_ex, 4: L2 divu, 5: L2 divu, 6: L2 sigma, 7: L2 sigma_ex

    if (!HasExactSol())
        DebugStop();

    errors.Resize(NEvalErrors());

    TPZManVector<STATE, 4> sol_exact(4);
    TPZFNMatrix<9, STATE> gradsol_exact(3, 3);

    // Getting the exact solution for velocity, pressure and velocity gradient
    fExactSol(data[EVindex].x, sol_exact, gradsol_exact);
    STATE p_exact = sol_exact[3];

    // Getting the numeric solution for velocity, pressure, and velocity gradient
    TPZManVector<STATE> u_h(3, 0.0);
    TPZManVector<STATE> p_h(1, 0.0);
    TPZFNMatrix<10, STATE> gradv_h = data[EVindex].dsol[0];

    this->Solution(data, VariableIndex("Velocity"), u_h);
    this->Solution(data, VariableIndex("Pressure"), p_h);

    STATE diff_v, diff_p, diff_div;

    // Pressure error
    diff_p = p_h[0] - p_exact;
    errors[0] = diff_p * diff_p;
    errors[1] = p_exact * p_exact;

    // Velocity Error
    errors[2] = 0.0;
    errors[3] = 0.0;
    for (int i = 0; i < fDimension; i++)
    {
        diff_v = u_h[i] - sol_exact[i];
        errors[2] += diff_v * diff_v;
        errors[3] += sol_exact[i] * sol_exact[i];
    }

    // Velocity Divergence Error
    STATE div_exact = 0.0, div_h = 0.0;
    for (int i = 0; i < fDimension; i++)
    {
        div_exact += gradsol_exact(i, i);
        div_h += gradv_h(i, i);
    }

    diff_div = div_h - div_exact;
    errors[4] = diff_div * diff_div;
    errors[5] = div_exact * div_exact;

    // Stress Tensor Error
    const int n = fDimension * (fDimension + 1) / 2;
    TPZFNMatrix<6, REAL> sigma_exact(n, 1), sigma_h(n, 1);

    StressTensor(gradsol_exact, sigma_exact, p_exact);
    StressTensor(gradv_h, sigma_h, p_h[0]);

    errors[6] = 0.0;
    errors[7] = 0.0;
    for (int i = 0; i < fDimension; i++)
    {
        const STATE diff_sigma = sigma_h(i, 0) - sigma_exact(i, 0);
        errors[6] += diff_sigma * diff_sigma;
        errors[7] += sigma_exact(i, 0) * sigma_exact(i, 0);
    }

    for (int i = fDimension; i < n; i++)
    {
        const STATE diff_sigma = sigma_h(i, 0) - sigma_exact(i, 0);
        errors[6] += 2.0 * diff_sigma * diff_sigma;
        errors[7] += 2.0 * sigma_exact(i, 0) * sigma_exact(i, 0);
    }

    // Deviatoric Stress Tensor
    TPZFNMatrix<6, REAL> devSigma_exact(n, 1), devSigma_h(n, 1);

    DeviatoricStressTensor(gradsol_exact, devSigma_exact);
    DeviatoricStressTensor(gradv_h, devSigma_h);

    errors[8] = 0.0;
    errors[9] = 0.0;
    for (int i = 0; i < fDimension; i++)
    {
        const STATE diff_devsigma = devSigma_h(i, 0) - devSigma_exact(i, 0);
        errors[8] += diff_devsigma * diff_devsigma;
        errors[9] += devSigma_exact(i, 0) * devSigma_exact(i, 0);
    }

    for (int i = fDimension; i < n; i++)
    {
        const STATE diff_devsigma = devSigma_h(i, 0) - devSigma_exact(i, 0);
        errors[8] += 2.0 * diff_devsigma * diff_devsigma;
        errors[9] += 2.0 * devSigma_exact(i, 0) * devSigma_exact(i, 0);
    }
}

void TPZHybridStokes::StressTensor(const TPZFNMatrix<10, STATE> &gradU, TPZFNMatrix<6, REAL> &sigma, REAL pressure)
{

    const int n = fDimension * (fDimension + 1) / 2;
    TPZFNMatrix<6, REAL> strain(n, 1, 0.0);
    TPZFNMatrix<36, REAL> D(n, n, 0.0);

    StrainTensor(gradU, strain);
    ViscosityTensor(D);

    D.Multiply(strain, sigma);

    for (int i = 0; i < fDimension; i++)
    {
        sigma(i, 0) -= pressure;
    }
}

void TPZHybridStokes::DeviatoricStressTensor(const TPZFNMatrix<10, STATE> &gradU, TPZFNMatrix<6, REAL> &sigma)
{

    const int n = fDimension * (fDimension + 1) / 2;
    TPZFNMatrix<6, REAL> strain(n, 1, 0.0);
    TPZFNMatrix<36, REAL> D(n, n, 0.0);

    StrainTensor(gradU, strain);
    ViscosityTensor(D);

    D.Multiply(strain, sigma);
}

void TPZHybridStokes::StrainTensor(const TPZFNMatrix<10, STATE> &gradU, TPZFNMatrix<6, REAL> &epsilon)
{

    const int n = fDimension * (fDimension + 1) / 2;
    int cont = fDimension - 1;

    for (int i = 0; i < fDimension; i++)
    {

        epsilon(i, 0) = gradU(i, i);

        for (int j = i + 1; j < fDimension; j++)
        {
            epsilon(++cont, 0) = 0.5 * (gradU(i, j) + gradU(j, i));
        }
    }
}

void TPZHybridStokes::ViscosityTensor(TPZFNMatrix<36, REAL> &D)
{
    int n = fDimension * (fDimension + 1) / 2;

    for (int i = 0; i < n; i++)
    {
        D(i, i) = 2 * fViscosity;
    }
}

void TPZHybridStokes::FromVoigt(const TPZFNMatrix<6, STATE> &sigmaVoigt, TPZFNMatrix<9, STATE> &sigma) const
{
    int cont = fDimension - 1;

    for (int i = 0; i < fDimension; i++)
    {
        sigma(i, i) = sigmaVoigt(i, 0);

        for (int j = i + 1; j < fDimension; j++)
        {
            sigma(i, j) = sigmaVoigt(++cont, 0);
            sigma(j, i) = sigmaVoigt(cont, 0);
        }
    }
}