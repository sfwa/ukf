#ifndef UKF_COVARIANCE_H
#define UKF_COVARIANCE_H

#include <Eigen/Cholesky>
#include <Eigen/Eigenvalues>
#include <Eigen/Householder>
#include <algorithm>
#include <limits>

namespace UKF {
namespace Detail {

/*
Return the lower triangular factor L, with a nonnegative diagonal, of the sum
of the outer products of the rows of A: L L^T = A^T A. Its transpose is the
triangular factor of the QR decomposition of A.

Each column of A in turn is reflected onto its diagonal entry by a
Householder reflection I - tau v v^T, where v is one followed by the
essential part which is left below the diagonal. The reflection is applied
to one column at a time, which for matrices as small as a filter's takes
about half as long as Eigen's HouseholderQR. Any non-finite entry of A
reaches L.
*/
template <typename Derived>
auto lower_factor(const Eigen::MatrixBase<Derived>& vectors) {
    using Scalar = typename Derived::Scalar;
    using Index = typename Derived::Index;
    Eigen::Matrix<Scalar, Derived::RowsAtCompileTime, Derived::ColsAtCompileTime, Eigen::ColMajor,
        Derived::MaxRowsAtCompileTime, Derived::MaxColsAtCompileTime> a = vectors;
    Eigen::Matrix<Scalar, Derived::ColsAtCompileTime, Derived::ColsAtCompileTime, Eigen::ColMajor,
        Derived::MaxColsAtCompileTime, Derived::MaxColsAtCompileTime> factor;
    factor.setZero(a.cols(), a.cols());
    for (Index k = 0; k < std::min(a.rows(), a.cols()); ++k) {
        const Index rows = a.rows() - k;
        Scalar tau = 0, pivot = 0;
        a.col(k).tail(rows).makeHouseholderInPlace(tau, pivot);
        a(k, k) = Scalar(1);
        for (Index j = k + 1; j < a.cols(); ++j) {
            a.col(j).tail(rows) -= (tau * a.col(k).tail(rows).dot(a.col(j).tail(rows))) * a.col(k).tail(rows);
        }

        /* Row k of the triangular factor is column k of L, up to its sign. */
        a(k, k) = pivot;
        factor.col(k).tail(a.cols() - k) = (pivot < Scalar(0) ? Scalar(-1) : Scalar(1)) *
            a.row(k).tail(a.cols() - k).transpose();
    }
    return factor;
}

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
    const MatrixType triangular = lower_factor(factor.transpose());
    if (!triangular.allFinite()) {
        return false;
    }
    root = triangular;
    return true;
}

}
}
#endif
