// This very slow. I stopped to computation after several minutes. Its definitely not zero cost abstraction.
#include "mex.hpp"
#include "mexAdapter.hpp"

using namespace matlab::data;
using matlab::mex::ArgumentList;

class MexFunction : public matlab::mex::Function {
    public:
    void operator()(ArgumentList outputs, ArgumentList inputs) {
        std::shared_ptr<matlab::engine::MATLABEngine> matlabPtr = getEngine();
        ArrayFactory factory;
        checkArguments(outputs, inputs);
        TypedArray<double> A = std::move(inputs[0]);
        TypedArray<double> b = std::move(inputs[1]);
        double tol = inputs[2][0];
        size_t n = A.getDimensions()[0];
        size_t nb = b.getDimensions()[1];
        size_t nnb = n + nb;
        size_t n1 = n + 1;
        TypedArray<double> Ab = factory.createArray<double>({n, nnb});
        for (size_t i = 0; i < n; i++) {
            for (size_t j = 0; j < n; j++) {
                Ab[i][j] = A[i][j];
            }
            for (size_t j = 0; j < nb; j++) {
                Ab[i][j+n] = b[i][j];
            }
        }
        for (size_t k = 0; k < n-1; k++) {
            size_t j = k;
            double maxVal = std::abs(Ab[k][k]);
            for (size_t i = k+1; i < n; i++) {
                double val = std::abs(Ab[i][k]);
                if (val > maxVal) {
                    j = i;
                    maxVal = val;
                }
            }
            if (j != k) {
                for (size_t i = 0; i < nnb; i++) {
                    double temp = Ab[k][i];
                    Ab[k][i] = Ab[j][i];
                    Ab[j][i] = temp;
                }
            }
            if (Ab[k][k] == 0) {
                matlabPtr->feval(u"error", 0, std::vector<Array>({factory.createScalar("Matrix is singular.")}));
            } else {
                for (size_t i = k+1; i < n; i++) {
                    double l = Ab[i][k] / Ab[k][k];
                    for (size_t j = k+1; j < nnb; j++) {
                        Ab[i][j] -= l * Ab[k][j];
                    }
                }
            }
        }
        double maxDiag = std::abs(Ab[0][0]);
        for (size_t i = 1; i < n; i++) {
            double val = std::abs(Ab[i][i]);
            if (val > maxDiag) {
                maxDiag = val;
            }
        }
        for (size_t i = 0; i < n; i++) {
            if (std::abs(Ab[i][i]) / maxDiag < tol) {
                matlabPtr->feval(u"warning", 0, std::vector<Array>({factory.createScalar("Matrix is close to singular or badly scaled. Results may be inaccurate.")}));
                break;
            }
        }
        for (int i = n-1; i >= 0; i--) {
            for (size_t j = 0; j < nb; j++) {
                Ab[i][j+n1] /= Ab[i][i];
            }
            for (int j = i-1; j >= 0; j--) {
                for (size_t k = 0; k < nb; k++) {
                    Ab[j][k+n1] -= Ab[j][i] * Ab[i][k+n1];
                }
            }
        }
        TypedArray<double> x = factory.createArray<double>({n, nb});
        for (size_t i = 0; i < n; i++) {
            for (size_t j = 0; j < nb; j++) {
                x[i][j] = Ab[i][j+n1];
            }
        }
        outputs[0] = x;
    }
    void checkArguments(ArgumentList outputs, ArgumentList inputs) {
        std::shared_ptr<matlab::engine::MATLABEngine> matlabPtr = getEngine();
        ArrayFactory factory;
        if (inputs.size() != 3) {
            matlabPtr->feval(u"error", 0, std::vector<Array>({factory.createScalar("Three input arguments required.")}));
        }
        if (inputs[0].getType() != ArrayType::DOUBLE || inputs[0].getType() == ArrayType::COMPLEX_DOUBLE) {
            matlabPtr->feval(u"error", 0, std::vector<Array>({factory.createScalar("Input matrix must be type double.")}));
        }
        if (inputs[1].getType() != ArrayType::DOUBLE || inputs[1].getType() == ArrayType::COMPLEX_DOUBLE) {
            matlabPtr->feval(u"error", 0, std::vector<Array>({factory.createScalar("Input vector must be type double.")}));
        }
        if (inputs[2].getType() != ArrayType::DOUBLE || inputs[2].getType() == ArrayType::COMPLEX_DOUBLE) {
            matlabPtr->feval(u"error", 0, std::vector<Array>({factory.createScalar("Input tolerance must be type double.")}));
        }
        if (outputs.size() > 1) {
            matlabPtr->feval(u"error", 0, std::vector<Array>({factory.createScalar("Too many output arguments.")}));
        }
    }
};
