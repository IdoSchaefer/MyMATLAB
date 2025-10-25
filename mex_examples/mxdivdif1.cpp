#include <cmath>
#include <iostream>
#include "mex.h"
/* Syntax: [polcoef, diagonal] = divdif(z, fz)
 *          polcoef = divdif(z, fz)
 *
 *
 */
using namespace std;
void compute_diag_real(double diagonal[], double z[], int coefi)
{
    int dtermi;
    for (dtermi = coefi - 1; dtermi>=0; dtermi--)
    {
        diagonal[dtermi] = (diagonal[dtermi + 1] - diagonal[dtermi])/(z[coefi] - z[dtermi]);
    }    
}

void compute_diag_realz(double Re_diagonal[], double Im_diagonal[], double z[], int coefi)
{
    int dtermi;
    double dz;
    for (dtermi = coefi - 1; dtermi>=0; dtermi--)
    {
        dz = z[coefi] - z[dtermi];
        Re_diagonal[dtermi] = (Re_diagonal[dtermi + 1] - Re_diagonal[dtermi])/dz;
        Im_diagonal[dtermi] = (Im_diagonal[dtermi + 1] - Im_diagonal[dtermi])/dz;
    }    
}

void compute_diag_comp(double Re_diagonal[], double Im_diagonal[], double Re_z[], double Im_z[], int coefi)
{
    int dtermi;
    double Re_dz, Im_dz, norm_dz, ezer;
    for (dtermi = coefi - 1; dtermi>=0; dtermi--)
    {
        Re_diagonal[dtermi] = (Re_diagonal[dtermi + 1] - Re_diagonal[dtermi]);
        Im_diagonal[dtermi] = (Im_diagonal[dtermi + 1] - Im_diagonal[dtermi]);
        Re_dz = Re_z[coefi] - Re_z[dtermi];
        Im_dz = Im_z[coefi] - Im_z[dtermi];
        norm_dz = pow(Re_dz, 2) + pow(Im_dz, 2);
        ezer = Re_diagonal[dtermi];
        Re_diagonal[dtermi] = (Re_diagonal[dtermi]*Re_dz + Im_diagonal[dtermi]*Im_dz)/norm_dz;
        Im_diagonal[dtermi] = (Im_diagonal[dtermi]*Re_dz - ezer*Im_dz)/norm_dz;                
    }    
}

/*double cmplx_division(Re_numerator, Im_numerator, Re_deno, Im_deno)
{
    
}*/
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    /*Macros for the output and input arguments: */
    #define mxpolcoef plhs[0]
    #define mxdiagonal plhs[1]
    #define mxz prhs[0]
    #define mxfz prhs[1]
    bool complex_z, complex_fz, complex_output;
    complex_z = mxIsComplex(mxz);
    complex_fz = mxIsComplex(mxfz);
    complex_output = (complex_z || complex_fz);
    double *Re_polcoef, *Re_diagonal, *Re_z, *Re_fz,
        *Im_polcoef, *Im_diagonal, *Im_z, *Im_fz;
    int dim, Npoints, dimi, coefi, dtermi;
    Re_z = mxGetPr(mxz);
    if (complex_z)
        Im_z = mxGetPi(mxz);
    Re_fz = mxGetPr(mxfz);
    if (complex_fz)
        Im_fz = mxGetPi(mxfz);
    dim = mxGetM(mxfz);
    Npoints = mxGetN(mxfz);
//    mxpolcoef = mxCreateDoubleMatrix(dim, Npoints, mxREAL);    
    if (complex_output)
        mxpolcoef = mxCreateDoubleMatrix(0, 0, mxCOMPLEX);
    else
        mxpolcoef = mxCreateDoubleMatrix(0, 0, mxREAL);    
    mxSetM(mxpolcoef, dim);
    mxSetN(mxpolcoef, Npoints);
    mxSetData(mxpolcoef, mxMalloc(sizeof(double)*dim*Npoints));
    Re_polcoef = mxGetPr(mxpolcoef);
    if (complex_output)
    {
        mxSetImagData(mxpolcoef, mxMalloc(sizeof(double)*dim*Npoints));
        Im_polcoef = mxGetPi(mxpolcoef);
    }
    if (nlhs < 2)
    {
        Re_diagonal = new double [Npoints];
        if (complex_output)
            Im_diagonal = new double [Npoints];
        for (dimi = 0; dimi<dim; dimi++)
        {
            Re_polcoef[dimi] = Re_fz[dimi];
            Re_diagonal[0] = Re_fz[dimi];
            if (complex_output)
            {
                Im_polcoef[dimi] = Im_fz[dimi];
                Im_diagonal[0] = Im_fz[dimi];               
            }
            for (coefi = 1; coefi<Npoints; coefi++)
            {
                Re_diagonal[coefi] = Re_fz[coefi*dim + dimi];
                if (complex_output)
                {
                    Im_diagonal[coefi] = Im_fz[coefi*dim + dimi];
                    if (complex_z)
                        compute_diag_comp(Re_diagonal, Im_diagonal, Re_z, Im_z, coefi);
                    else
                        compute_diag_realz(Re_diagonal, Im_diagonal, Re_z, coefi);
                    Im_polcoef[coefi*dim + dimi] = Im_diagonal[0];
                }
                else
                    compute_diag_real(Re_diagonal, Re_z, coefi);
                Re_polcoef[coefi*dim + dimi] = Re_diagonal[0];
            }
        }
        delete [] Re_diagonal;
        if (complex_output)
            delete [] Im_diagonal;
    }
    else
    // !!!To be completed later!!!
    {
//        mxdiagonal = mxCreateDoubleMatrix(dim, Npoints, mxREAL);
        mxdiagonal = mxCreateDoubleMatrix(dim, Npoints, mxREAL);
        mxSetM(mxdiagonal, dim);
        mxSetN(mxdiagonal, Npoints);
        mxSetData(mxdiagonal, mxMalloc(sizeof(double)*dim*Npoints));
        Re_diagonal = mxGetPr(mxdiagonal);        
        for (dimi = 0; dimi<dim; dimi++)
        {
            Re_polcoef[dimi] = Re_fz[dimi];
            Re_diagonal[dimi] = Re_fz[dimi];
            for (coefi = 1; coefi<Npoints; coefi++)
            {
                Re_diagonal[coefi*dim + dimi] = Re_fz[coefi*dim + dimi];
                for (dtermi = coefi - 1; dtermi>=0; dtermi--)
                {
                    Re_diagonal[dtermi*dim + dimi] = (Re_diagonal[(dtermi + 1)*dim + dimi] - Re_diagonal[dtermi*dim + dimi])/(Re_z[coefi] - Re_z[dtermi]);
                }
                Re_polcoef[coefi*dim + dimi] = Re_diagonal[dimi];
            }           
        }
    }
}