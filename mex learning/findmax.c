#include "mex.h"
#include "blas.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
  // Setup inputs
  const mxDouble *x = mxGetDoubles(prhs[0]);
  const mwIndex column = (mwIndex) mxGetScalar(prhs[1]) - 1;
  const mwSize m = mxGetM(prhs[0]);
  const mwIndex increment = 1;

  // Setup outputs
  plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
  mxDouble *index = mxGetDoubles(plhs[0]);
  
  // Call the Fortran BLAS function idamax_
  *index = (mxDouble)idamax_(&m, x + column*m, &increment);
}