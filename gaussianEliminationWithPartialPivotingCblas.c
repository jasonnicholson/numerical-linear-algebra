// Must be compiled with "mex -R2018a -lmwblas ..." 

#include "mex.h"
#include <math.h>
#include "blas.h"
#include <string.h>

// Detect the OS
#if defined(_WIN32) || defined(_WIN64)
    #define WINDOWS_OS
#elif defined(__unix__) || defined(__unix) || defined(__APPLE__)
    #define UNIX_OS
#endif

// Override function names for Windows
#ifdef WINDOWS_OS
    #define idamax_ idamax
    #define dswap_ dswap
    #define dscal_ dscal
    #define dger_ dger
    #define dtrsm_ dtrsm
#endif

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

    // Combining A and b. This aggregates the data needed in cache for the algorithm.
    mxArray *Ab_mxArray = mxCreateDoubleMatrix(nRowsA, nColumnsAandb, mxREAL);
    mxDouble *Ab = mxGetDoubles(Ab_mxArray);
    memcpy(Ab, A, nRowsA*nRowsA*sizeof(mxDouble));
    memcpy(Ab + nRowsA*nRowsA, b, nRowsA*nColumnsb*sizeof(mxDouble));
    
    // Perform gaussian elimination with partial pivoting
    mwSignedIndex ROW_INCREMENT = 1; // defines how big of step to take through the column vector.
    mwSignedIndex COLUMN_INCREMENT = nRowsA; // defines how big of step to take through the row vector.
    mwSignedIndex leadingDimensionOfAb = nRowsA; // defines the leading dimension of the matrix Ab.
    for (mwIndex kColumn = 0; kColumn < nColumnsA-1; kColumn++) {
        mwIndex kRow = kColumn;
        //Find the pivot index by looking down the column
        // Call BLAS function idamax_ . idamax_ returns 1-based indexing.
        mwSignedIndex nRowsInKcolumn = nRowsA - (mwSignedIndex)kRow;
        mwSignedIndex jRow = idamax_(&nRowsInKcolumn, Ab + kColumn*nRowsA + kRow, &ROW_INCREMENT) - 1 + kRow; // Subtract 1 to get 0-based indexing

        // Swap the k-th and j-th rows in Ab
        // Use dswap_ routine from BLAS
        dswap_(&nColumnsAandb, Ab + jRow, &nRowsA, Ab + kRow, &nRowsA);

        // Perform the gaussian elimination
        if (Ab[kRow + kColumn*nRowsA] == 0) {
            mxDestroyArray(Ab_mxArray);
            mexErrMsgIdAndTxt("MATLAB:gaussianEliminationWithPartialPivoting:matrixIsSingular", "Matrix is singular.");
        } else {
            // Calculate lower triangular part. This calculates the column of L matrix at the kth column.
            // Use dscal_ routine from BLAS
            mwSignedIndex nRowsInKcolumnMinus1 = nRowsInKcolumn - 1;
            mxDouble pivotInverse = 1.0/Ab[kRow + kColumn*nRowsA];
            mxDouble *lVector = Ab + (kRow + 1) + kColumn*nRowsA;
            dscal_(&nRowsInKcolumnMinus1, &pivotInverse, lVector, &ROW_INCREMENT);

            // Calculate upper triangular part. This calculates the entire submatrix of the augemented U matrix at the current kth step.
            // This handles Ly = b part of the algorithm as well.
            // Use dger_ routine from BLAS
            mwSignedIndex nRowsInSubMatrix = nRowsInKcolumnMinus1;
            mwSignedIndex nColumnsInSubMatrix = nColumnsAandb - (kColumn + 1);
            mxDouble alpha = -1.0; // defines the scalar alpha in the dger_ routine.
            mxDouble *uRowVector = Ab + kRow  + (kColumn+1) * nRowsA;
            mxDouble *uSubMatrix = Ab + (kRow + 1) + (kColumn+1) * nRowsA;
            dger_(&nRowsInSubMatrix, &nColumnsInSubMatrix, &alpha, lVector, &ROW_INCREMENT, uRowVector, &COLUMN_INCREMENT, uSubMatrix, &leadingDimensionOfAb);
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
    // use dtrsm_ routine from BLAS
    char side = 'L', uplo = 'U', transa = 'N', diag = 'N';
    mxDouble alpha = 1.0;
    mwSignedIndex leadingDimensionOfb = nRowsA;
    dtrsm_(&side, &uplo, &transa, &diag, &nRowsA, &nColumnsb, &alpha, Ab, &leadingDimensionOfAb, Ab + nRowsA * nColumnsA, &leadingDimensionOfb);

    // Copy the solution to the output
    plhs[0] = mxCreateDoubleMatrix(nRowsA, nColumnsb, mxREAL);
    mxDouble *x = mxGetDoubles(plhs[0]);
    memcpy(x, Ab + nRowsA * nColumnsA, nRowsA * nColumnsb * sizeof(mxDouble));

    // Clean up
    mxDestroyArray(Ab_mxArray);
}
