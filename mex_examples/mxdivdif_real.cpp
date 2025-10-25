#include <cmath>
#include <iostream>
#include "mex.h"
/* Syntax: [polcoef, diagonal] = divdif(z, fz)
 *          polcoef = divdif(z, fz)
 *
 *
 */
using namespace std;
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    /*Macros for the output and input arguments: */
    #define mxpolcoef plhs[0]
    #define mxdiagonal plhs[1]
    #define mxz prhs[0]
    #define mxfz prhs[1]
    double *polcoef, *diagonal, *z, *fz;
    int dim, Npoints, dimi, coefi, dtermi;
    z = mxGetPr(mxz);
    fz = mxGetPr(mxfz);
    dim = mxGetM(mxfz);
    Npoints = mxGetN(mxfz);
    //mxpolcoef = mxCreateDoubleMatrix(dim, Npoints, mxREAL);
    mxpolcoef = mxCreateDoubleMatrix(0, 0, mxREAL);
    mxSetM(mxpolcoef, dim);
    mxSetN(mxpolcoef, Npoints);
    mxSetData(mxpolcoef, mxMalloc(sizeof(double)*dim*Npoints));
    polcoef = mxGetPr(mxpolcoef);
    if (nlhs < 2)
    {
        diagonal = new double [Npoints];
        for (dimi = 0; dimi<dim; dimi++)
        {
            polcoef[dimi] = fz[dimi];
            diagonal[0] = fz[dimi];
            for (coefi = 1; coefi<Npoints; coefi++)
            {
                diagonal[coefi] = fz[coefi*dim + dimi];
                for (dtermi = coefi - 1; dtermi>=0; dtermi--)
                {
                    diagonal[dtermi] = (diagonal[dtermi + 1] - diagonal[dtermi])/(z[coefi] - z[dtermi]);
                }
                polcoef[coefi*dim + dimi] = diagonal[0];
            }
        }
        delete [] diagonal;
    }
    else
    {
        //mxdiagonal = mxCreateDoubleMatrix(dim, Npoints, mxREAL);
        mxdiagonal = mxCreateDoubleMatrix(dim, Npoints, mxREAL);
        mxSetM(mxdiagonal, dim);
        mxSetN(mxdiagonal, Npoints);
        mxSetData(mxdiagonal, mxMalloc(sizeof(double)*dim*Npoints));
        diagonal = mxGetPr(mxdiagonal);        
        for (dimi = 0; dimi<dim; dimi++)
        {
            polcoef[dimi] = fz[dimi];
            diagonal[dimi] = fz[dimi];
            for (coefi = 1; coefi<Npoints; coefi++)
            {
                diagonal[coefi*dim + dimi] = fz[coefi*dim + dimi];
                for (dtermi = coefi - 1; dtermi>=0; dtermi--)
                {
                    diagonal[dtermi*dim + dimi] = (diagonal[(dtermi + 1)*dim + dimi] - diagonal[dtermi*dim + dimi])/(z[coefi] - z[dtermi]);
                }
                polcoef[coefi*dim + dimi] = diagonal[dimi];
            }           
        }
    }
}