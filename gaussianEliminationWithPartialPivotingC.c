#include "mex.h"
#include <math.h>

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    double *A, *b, *x, *Ab, tol;
    size_t n, nb, nnb;
    mwIndex i, j, k;
    mxArray *Ab_mxArray;

    /* Check for proper number of input and output arguments */
    if (nrhs != 3) {
        mexErrMsgIdAndTxt("MATLAB:gaussianEliminationWithPartialPivoting:invalidNumInputs", "Three input arguments required.");
    }
    if (nlhs > 1) {
        mexErrMsgIdAndTxt("MATLAB:gaussianEliminationWithPartialPivoting:maxlhs", "Too many output arguments.");
    }

    /* Check input argument types */
    if (!mxIsDouble(prhs[0]) || 
        !mxIsDouble(prhs[1]) || 
        !mxIsDouble(prhs[2])) {
        mexErrMsgIdAndTxt("MATLAB:gaussianEliminationWithPartialPivoting:inputNotRealScalarDouble", "Input arguments must be real double.");
    }

    /* Get the size of the matrix A and b vector */
    n = mxGetM(prhs[0]);
    nb = mxGetN(prhs[1]);
    nnb = n + nb;

    /* Combining A and b signifies most of the data needed in cache for processing */
    Ab_mxArray = mxCreateDoubleMatrix(n, nnb, mxREAL);
    Ab = mxGetDoubles(Ab_mxArray);
    A = mxGetDoubles(prhs[0]);
    b = mxGetDoubles(prhs[1]);
    tol = mxGetScalar(prhs[2]);
    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            Ab[i + j*n] = A[i + j*n];
        }
        for (j = 0; j < nb; j++) {
            Ab[i + (j+n)*n] = b[i + j*n];
        }
    }

    /* Perform gaussian elimination with partial pivoting */
    for (k = 0; k < n-1; k++) {
        /* Find the pivot index */
        j = k;
        for (i = k+1; i < n; i++) {
            if (fabs(Ab[i + k*n]) > fabs(Ab[j + k*n])) {
                j = i;
            }
        }

        /* Swap the k-th and j-th rows in Ab */
        for (i = 0; i < nnb; i++) {
            double temp = Ab[k + i*n];
            Ab[k + i*n] = Ab[j + i*n];
            Ab[j + i*n] = temp;
        }

        /* Perform the gaussian elimination */
        if (Ab[k + k*n] == 0) {
            mexErrMsgIdAndTxt("MATLAB:gaussianEliminationWithPartialPivoting:matrixIsSingular",
                              "Matrix is singular.");
        } else {
            for (i = k+1; i < n; i++) {
                double lVector = Ab[i + k*n] / Ab[k + k*n]; /* Calculate lower triangular part */
                for (j = k+1; j < nnb; j++) {
                    Ab[i + j*n] -= lVector * Ab[k + j*n]; /* Calculate upper triangular part */
                }
            }
        }
    }

    /* Check for rank deficiency or near singularity */
    mxDouble maxDiag = 0.0;
    for (i = 0; i < n; i++) {
      mxDouble currentDiag = fabs(Ab[i + i*n]);
      if (currentDiag > maxDiag) {
        maxDiag = currentDiag;
      }
    }
    for (i = 0; i < n; i++) {
        if (fabs(Ab[i + i*n])/maxDiag < tol) {
            mexWarnMsgIdAndTxt("MATLAB:gaussianEliminationWithPartialPivoting:matrixIsCloseToSingular", "Matrix is close to singular or badly scaled. Results may be inaccurate.");
            break;
        }
    }

    /* Backward substitution to solve Ux = y */
    plhs[0] = mxCreateDoubleMatrix(n, nb, mxREAL);
    x = mxGetPr(plhs[0]);
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

    /* Clean up */
    mxDestroyArray(Ab_mxArray);
}
