// Must be compiled with "mex -R2018a -lmwblas ..." 

#include "mex.h"
#include <math.h>
#include "blas.h"
#include <string.h>

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    mxDouble *A, *b, *x, *Ab, *tol;
    mwSize nb, nnb;
    mwIndex i, j, k;
    mwSignedIndex n;
    mxArray *Ab_mxArray;

  // Get the size of the matrix A and b vector
    n = mxGetM(prhs[0]);
    nb = mxGetN(prhs[1]);
    nnb = n + nb;

    // Check for proper number of input and output arguments
    if (nrhs != 3) {
        mexErrMsgIdAndTxt("MATLAB:gaussianEliminationWithPartialPivoting:invalidNumInputs",
                          "Three input arguments required.");
    }
    if (nlhs > 1) {
        mexErrMsgIdAndTxt("MATLAB:gaussianEliminationWithPartialPivoting:maxlhs",
                          "Too many output arguments.");
    }
    // Check input argument types
    if (!mxIsDouble(prhs[0]) || 
        !mxIsDouble(prhs[1]) || 
        !mxIsDouble(prhs[2])) {
        mexErrMsgIdAndTxt("MATLAB:gaussianEliminationWithPartialPivoting:inputNotRealScalarDouble",
                          "Input arguments must be real double.");
    }
    // Check that A and b are compatible
    if (mxGetN(prhs[0]) != n) {
        mexErrMsgIdAndTxt("MATLAB:gaussianEliminationWithPartialPivoting:invalidInput",
                          "Matrix A must be square.");
    }
    if (mxGetM(prhs[1]) != n) {
        mexErrMsgIdAndTxt("MATLAB:gaussianEliminationWithPartialPivoting:invalidInput",
                          "Matrix A and vector b must have the same number of rows.");
    }

    // Get the inputs
    A = mxGetDoubles(prhs[0]);
    b = mxGetDoubles(prhs[1]);
    tol = mxGetDoubles(prhs[2]);

    // Combining A and b signifies most of the data needed in cache for processing
    Ab_mxArray = mxCreateDoubleMatrix(n, nnb, mxREAL);
    Ab = mxGetDoubles(Ab_mxArray);
    memcpy(Ab, A, n*n*sizeof(mxDouble));
    memcpy(Ab + n*n, b, n*nb*sizeof(mxDouble));
    
    // Perform gaussian elimination with partial pivoting
    for (k = 0; k < n-1; k++) {
        //Find the pivot index by looking down the column
        ptrdiff_t increment = 1; // defines how big of step to take through the vector.
        // Call BLAS function idamax_ . idamax_ returns 1-based indexing.
        ptrdiff_t nn = (ptrdiff_t)(n-k);
        j = idamax_(&nn, Ab + k*n + k, &increment) - 1 + k; // Subtract 1 to get 0-based indexing

        // Swap the k-th and j-th rows in Ab
        // Use dswap_ routine from BLAS
        dswap_(&nnb, Ab + j, &n, Ab + k, &n);

        // Perform the gaussian elimination
        if (Ab[k + k*n] == 0) {
            mexErrMsgIdAndTxt("MATLAB:gaussianEliminationWithPartialPivoting:matrixIsSingular", "Matrix is singular.");
        } else {
            for (i = k+1; i < n; i++) {
                double lVector = Ab[i + k*n] / Ab[k + k*n]; // Calculate lower triangular part
                for (j = k+1; j < nnb; j++) {
                    Ab[i + j*n] -= lVector * Ab[k + j*n]; // Calculate upper triangular part
                }
            }
        }
    }

    // Check for rank deficiency or near singularity
    mwSignedIndex maxDiagIndex = idamax_(&n, Ab, &n) - 1; // find the maximum diagonal element
    mxDouble maxDiag = fabs(Ab[maxDiagIndex]);
    for (i = 0; i < n; i++) {
        if (fabs(Ab[i + i*n])/maxDiag < *tol) {
            mexWarnMsgIdAndTxt("MATLAB:gaussianEliminationWithPartialPivoting:matrixIsCloseToSingular", "Matrix is close to singular or badly scaled. Results may be inaccurate.");
            break;
        }
    }

    // Backward substitution to solve Ux = y
    plhs[0] = mxCreateDoubleMatrix(n, nb, mxREAL);
    x = mxGetDoubles(plhs[0]);
    for (i = n; i-- > 0;) {
        for (j = 0; j < nb; j++) {
            x[i + j*n] = Ab[i + (j+n)*n] / Ab[i + i*n];
        }
        for (j = 0; j < i; j++) {
            for (k = 0; k < nb; k++) {
                Ab[j + (k+n)*n] -= Ab[j + i*n] * x[i + k*n];
            }
        }
    }

    // Clean up
    mxDestroyArray(Ab_mxArray);
}
