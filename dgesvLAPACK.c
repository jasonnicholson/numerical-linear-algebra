// Must be compiled with "mex -R2018a -lmwblas ..." 

#include "mex.h"
#include <math.h>
#include "lapack.h"
#include <string.h>

// Detect the OS
#if defined(_WIN32) || defined(_WIN64)
    #define WINDOWS_OS
#elif defined(__unix__) || defined(__unix) || defined(__APPLE__)
    #define UNIX_OS
#endif

// Override function names for Windows
#ifdef WINDOWS_OS
    #define dgesv_ dgesv
#endif

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
  // Get the size of the matrix A and b vector
  // Normally we would use mwSize, but we are using mwSignedIndex to avoid warnings with BLAS functions when using -Wall.
    mwSignedIndex nRowsA = mxGetM(prhs[0]);
    mwSignedIndex nColumnsb = mxGetN(prhs[1]);

    // Check for proper number of input and output arguments
    if (nrhs < 2 && nrhs > 3) {
        mexErrMsgIdAndTxt("MATLAB:dgesvLAPACK:invalidNumInputs", "2 or 3 input arguments required.");
    }
    if (nlhs > 1) {
        mexErrMsgIdAndTxt("MATLAB:dgesvLAPACK:maxlhs", "Too many output arguments.");
    }
    // Check input argument types
    if (!mxIsDouble(prhs[0]) || 
        !mxIsDouble(prhs[1])) {
        mexErrMsgIdAndTxt("MATLAB:dgesvLAPACK:inputNotRealDouble", "Input arguments must be real double.");
    }
    // Check that A and b are compatible
    if (mxGetN(prhs[0]) != nRowsA) {
        mexErrMsgIdAndTxt("MATLAB:dgesvLAPACK:invalidInput", "Matrix A must be square.");
    }
    if (mxGetM(prhs[1]) != nRowsA) {
        mexErrMsgIdAndTxt("MATLAB:dgesvLAPACK:invalidInput",  "Matrix A and vector b must have the same number of rows.");
    }

    // Copy A and b
    mxArray *A = mxDuplicateArray(prhs[0]);
    mxArray *b = mxDuplicateArray(prhs[1]);
    
    // needed for LAPACK
    mwSignedIndex *pivotIndex = (mwSignedIndex *)mxCalloc(nRowsA, sizeof(mwSignedIndex));
    if (pivotIndex == NULL) {
      mxFree(pivotIndex);
      mxDestroyArray(A);
      mexErrMsgTxt("Memory allocation for pivot indices failed.");
    }
    mwSignedIndex info; // tells of convergence

    // Call LAPACK dgesv
    dgesv_(&nRowsA, &nColumnsb, mxGetDoubles(A), &nRowsA, pivotIndex, mxGetDoubles(b), &nRowsA, &info);

    // Clean up
    mxFree(pivotIndex);
    mxDestroyArray(A);

    // Check the value of info and handle any errors
    if (info == 0) {
      // Copy the solution to the output
      plhs[0] = b;
    } else if (info < 0) {
        mexErrMsgIdAndTxt("MATLAB:dgesvLAPACK:InvalidArgument", "Argument %d to dgesv_ has an illegal value.", -info);
    } else if (info > 0) {
        mexErrMsgIdAndTxt("MATLAB:dgesvLAPACK:SingularMatrix", "Matrix A is singular. Solution could not be computed.");
    }
}
