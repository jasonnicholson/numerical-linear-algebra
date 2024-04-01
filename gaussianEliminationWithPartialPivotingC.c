// Must be compiled with "mex -R2018a -lmwblas ..." 

#include "mex.h"
#include <math.h>
#include "blas.h"
#include <string.h>

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
  // Get the size of the matrix A and b vector
  // Normally we would use mwSize, but we are using mwSignedIndex to avoid warnings with BLAS functions when using -Wall.
    mwSignedIndex nRowsA = mxGetM(prhs[0]);
    mwSignedIndex nColumnsA = nRowsA;
    mwSignedIndex nColumnsb = mxGetN(prhs[1]);
    mwSignedIndex nColumnsAandb = nRowsA + nColumnsb;

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
        mexErrMsgIdAndTxt("MATLAB:gaussianEliminationWithPartialPivoting:inputNotRealScalarDouble", "Input arguments must be real double.");
    }
    // Check that A and b are compatible
    if (mxGetN(prhs[0]) != nRowsA) {
        mexErrMsgIdAndTxt("MATLAB:gaussianEliminationWithPartialPivoting:invalidInput", "Matrix A must be square.");
    }
    if (mxGetM(prhs[1]) != nRowsA) {
        mexErrMsgIdAndTxt("MATLAB:gaussianEliminationWithPartialPivoting:invalidInput",  "Matrix A and vector b must have the same number of rows.");
    }

    // Get the inputs
    mxDouble *A = mxGetDoubles(prhs[0]);
    mxDouble *b = mxGetDoubles(prhs[1]);
    mxDouble tol = mxGetScalar(prhs[2]);

    // Combining A and b. This the main data needed in cache for the algorithm.
    mxArray *Ab_mxArray = mxCreateDoubleMatrix(nRowsA, nColumnsAandb, mxREAL);
    mxDouble *Ab = mxGetDoubles(Ab_mxArray);
    memcpy(Ab, A, nRowsA*nRowsA*sizeof(mxDouble));
    memcpy(Ab + nRowsA*nRowsA, b, nRowsA*nColumnsb*sizeof(mxDouble));
    
    // Perform gaussian elimination with partial pivoting
    for (mwIndex kColumn = 0; kColumn < nColumnsA-1; kColumn++) {
        mwIndex kRow = kColumn;
        //Find the pivot index by looking down the column
        mwSignedIndex increment = 1; // defines how big of step to take through the vector.
        // Call BLAS function idamax_ . idamax_ returns 1-based indexing.
        mwSignedIndex nRowsInKcolumn = nRowsA - (mwSignedIndex)kRow;
        mwSignedIndex jRow = idamax_(&nRowsInKcolumn, Ab + kColumn*nRowsA + kRow, &increment) - 1 + kRow; // Subtract 1 to get 0-based indexing

        // Swap the k-th and j-th rows in Ab
        // Use dswap_ routine from BLAS
        dswap_(&nColumnsAandb, Ab + jRow, &nRowsA, Ab + kColumn, &nRowsA);

        // Perform the gaussian elimination
        if (Ab[kRow + kColumn*nRowsA] == 0) {
            mexErrMsgIdAndTxt("MATLAB:gaussianEliminationWithPartialPivoting:matrixIsSingular", "Matrix is singular.");
        } else {
            for (mwIndex iRow = kRow+1; iRow < nRowsA; iRow++) {
                Ab[iRow + kColumn*nRowsA] = Ab[iRow + kColumn*nRowsA] / Ab[kRow + kColumn*nRowsA]; // Calculate lower triangular part
                for (mwIndex jColumn = kColumn+1; jColumn < nColumnsAandb; jColumn++) {
                    Ab[iRow + jColumn*nRowsA] -= Ab[iRow + kColumn*nRowsA] * Ab[kRow + jColumn*nRowsA]; // Calculate upper triangular part
                }
            }
        }
    }

    // Check for rank deficiency or near singularity
    mwSignedIndex incrementBetweenDiagonalElements = nRowsA + 1; // defines how big of step to take through the vector.
    mwSignedIndex maxDiagIndex = idamax_(&nRowsA, Ab, &incrementBetweenDiagonalElements) - 1; // find the maximum diagonal element
    mxDouble maxDiag = fabs(Ab[maxDiagIndex + maxDiagIndex*nRowsA]);
    for (mwIndex iColumn = 0; iColumn < nColumnsA; iColumn++) {
        mwIndex iRow = iColumn;
        if (fabs(Ab[iRow + iColumn*nRowsA])/maxDiag < tol) {
            mexWarnMsgIdAndTxt("MATLAB:gaussianEliminationWithPartialPivoting:matrixIsCloseToSingular", "Matrix is close to singular or badly scaled. Results may be inaccurate.");
            break;
        }
    }

    // Backward substitution to solve Ux = y
    plhs[0] = mxCreateDoubleMatrix(nRowsA, nColumnsb, mxREAL);
    mxDouble *x = mxGetDoubles(plhs[0]);
    for (mwIndex iRow = nRowsA; iRow-- > 0;) {
        mwIndex iColumn = iRow;
        for (mwIndex jColumnb = 0; jColumnb < nColumnsb; jColumnb++) {
            x[iRow + jColumnb*nRowsA] = Ab[iRow + (jColumnb+nColumnsA)*nRowsA] / Ab[iRow + iColumn*nRowsA];
        }
        for (mwIndex jRow = 0; jRow < iRow; jRow++) {
            for (mwIndex kColumnb = 0; kColumnb < nColumnsb; kColumnb++) {
                Ab[jRow + (kColumnb+nColumnsA)*nRowsA] -= Ab[jRow + iColumn*nRowsA] * x[iRow + kColumnb*nRowsA];
            }
        }
    }

    // Clean up
    mxDestroyArray(Ab_mxArray);
}
