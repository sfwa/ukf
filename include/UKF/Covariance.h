#ifndef UKF_COVARIANCE_H
#define UKF_COVARIANCE_H

#include <Eigen/Cholesky>
#include <Eigen/Eigenvalues>
#include <Eigen/QR>
#include <algorithm>
#include <limits>
#include <stdexcept>

namespace UKF {
namespace Detail {

/*
Return a lower triangular factor of a positive semidefinite covariance.
Keep LLT for the usual positive definite case. For a singular covariance,
clip only roundoff-sized negative eigenvalues and triangularize its spectral
factor. Fixed-size state matrices keep all Eigen storage on the stack.
The optional reference scale accounts for cancellation in covariance updates.
*/
template <typename MatrixType>
MatrixType covariance_root(const MatrixType& covariance,
        typename MatrixType::Scalar reference_scale = 0) {
    using Scalar = typename MatrixType::Scalar;
    if (!covariance.allFinite()) {
        throw std::domain_error("Covariance matrix is not finite");
    }
    const Eigen::LLT<MatrixType> llt(covariance);
    if (llt.info() == Eigen::Success) {
        return llt.matrixL();
    }

    const Scalar tolerance = Scalar(16) * std::numeric_limits<Scalar>::epsilon() *
        Scalar(covariance.rows()) * std::max(reference_scale, covariance.cwiseAbs().maxCoeff());
    const Eigen::SelfAdjointEigenSolver<MatrixType> eig(covariance);
    if (eig.info() != Eigen::Success || eig.eigenvalues().minCoeff() < -tolerance) {
        throw std::domain_error("Covariance matrix is not positive semidefinite");
    }
    const MatrixType factor = eig.eigenvectors() * eig.eigenvalues().cwiseMax(Scalar(0)).cwiseSqrt().asDiagonal();
    MatrixType root = factor.transpose().householderQr().matrixQR().template triangularView<Eigen::Upper>();
    root.transposeInPlace();
    return root;
}

}
}
#endif
