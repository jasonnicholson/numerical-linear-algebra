#include "mex.h"
#include "blas.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
  // Setup inputs
  mxDouble *x = mxGetDoubles(prhs[0]);
  mwSize m = mxGetM(prhs[0]);
  mwIndex i = 1;

  // Setup outputs
  plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
  mxDouble *index = mxGetDoubles(plhs[0]);
  

  // Call the Fortran BLAS function idamax_
  *index = (mxDouble)idamax_(&m, x, &i);
}