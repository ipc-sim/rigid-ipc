#include "eigen_ext.hpp"

#include <Eigen/Eigenvalues>

#include <logger.hpp>

namespace Eigen {

// Matrix Projection onto Positive Semi-Definite Cone
MatrixXd project_to_psd(const MatrixXd& A)
{
    // https://math.stackexchange.com/q/2776803
    SelfAdjointEigenSolver<MatrixXd> eigensolver(A);
    if (eigensolver.info() != Success) {
        spdlog::error(
            "unable to project matrix onto positive semi-definite cone");
    }
    // Check if all eigen values are zero or positive.
    // The eigenvalues are sorted in increasing order.
    if (eigensolver.eigenvalues()[0] >= 0.0) {
        return A;
    }
    DiagonalMatrix<double, Dynamic> D(eigensolver.eigenvalues());
    // Save a little time and only project the negative values
    for (int i = 0; i < A.rows(); i++) {
        if (D.diagonal()[i] < 0.0) {
            D.diagonal()[i] = 0.0;
        } else {
            break;
        }
    }
    return eigensolver.eigenvectors() * D
        * eigensolver.eigenvectors().transpose();
}

} // namespace Eigen
