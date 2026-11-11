#ifndef UKF_COVARIANCE_H
#define UKF_COVARIANCE_H

#include <Eigen/Cholesky>
#include <Eigen/Eigenvalues>
#include <Eigen/QR>
#include <algorithm>
#include <limits>

namespace UKF {
namespace Detail {

/*
Compute a lower triangular factor of a positive semidefinite covariance.
Keep LLT for the usual positive definite case. For a singular covariance,
clip only roundoff-sized negative eigenvalues and triangularize its spectral
factor. Fixed-size state matrices keep all Eigen storage on the stack.
The optional reference scale accounts for cancellation in covariance updates.

Return false, leaving root unchanged, if the covariance is not finite or is
indefinite beyond roundoff. The library reports failures through return
values and never throws, so it can be built with exceptions disabled.
*/
template <typename MatrixType>
bool covariance_root(const MatrixType& covariance, MatrixType& root,
        typename MatrixType::Scalar reference_scale = 0) {
    using Scalar = typename MatrixType::Scalar;
    if (!covariance.allFinite()) {
        return false;
    }
    const Eigen::LLT<MatrixType> llt(covariance);
    if (llt.info() == Eigen::Success) {
        root = llt.matrixL();
        return true;
    }

    const Scalar tolerance = Scalar(16) * std::numeric_limits<Scalar>::epsilon() *
        Scalar(covariance.rows()) * std::max(reference_scale, covariance.cwiseAbs().maxCoeff());
    const Eigen::SelfAdjointEigenSolver<MatrixType> eig(covariance);
    if (eig.info() != Eigen::Success || eig.eigenvalues().minCoeff() < -tolerance) {
        return false;
    }
    const MatrixType factor = eig.eigenvectors() * eig.eigenvalues().cwiseMax(Scalar(0)).cwiseSqrt().asDiagonal();
    MatrixType triangular = factor.transpose().householderQr().matrixQR().template triangularView<Eigen::Upper>();
    triangular.transposeInPlace();
    if (!triangular.allFinite()) {
        return false;
    }
    root = triangular;
    return true;
}

}
}
#endif
