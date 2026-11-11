/* SPDX-License-Identifier: MIT */
#ifndef UKF_RANK_UPDATE_H
#define UKF_RANK_UPDATE_H

#include <Eigen/Core>
#include <cmath>
#include <type_traits>

namespace UKF {
namespace Detail {

/*
Update L L^T by weight * v v^T using orthogonal (positive weight) or
hyperbolic (negative weight) rotations. Keep the workspace bounded, unlike
Eigen's LLT rank-update temporary in unpatched Eigen releases. Only the lower
triangle is read or changed. Return -1 on success, or the index of the first
pivot a downdate would make nonpositive; callers must discard a partially
completed failed downdate. Zero pivots with zero update entries, as in
semidefinite factors, pass through unchanged.
*/
template <typename MatrixType, typename VectorType>
typename MatrixType::Index rank_update_lower(MatrixType& factor,
        const VectorType& vector, typename MatrixType::Scalar weight) {
    using Scalar = typename MatrixType::Scalar;
    using Index = typename MatrixType::Index;
    static_assert(std::is_floating_point<Scalar>::value,
        "UKF rank updates require real floating-point factors");
    static_assert(MatrixType::MaxRowsAtCompileTime != Eigen::Dynamic,
        "UKF rank updates require fixed or bounded matrix storage");
    eigen_assert(factor.rows() == factor.cols() && factor.rows() == vector.size());
    if (weight == Scalar(0)) {
        return Index(-1);
    }

    Eigen::Matrix<Scalar, MatrixType::RowsAtCompileTime, 1, 0,
        MatrixType::MaxRowsAtCompileTime, 1> work(factor.rows());
    work = std::sqrt(std::abs(weight)) * vector;

    for (Index k = 0; k < factor.rows(); ++k) {
        const Scalar diagonal = factor(k, k);
        const Scalar entry = work(k);
        if (weight > Scalar(0)) {
            const Scalar magnitude = std::hypot(diagonal, entry);
            if (magnitude == Scalar(0)) {
                continue;
            }
            // QR may supply negative diagonals; return a canonical positive pivot.
            const Scalar updated = magnitude;
            const Scalar cosine = diagonal / updated;
            const Scalar sine = entry / updated;
            factor(k, k) = updated;
            for (Index j = k + 1; j < factor.rows(); ++j) {
                const Scalar old_factor = factor(j, k);
                const Scalar old_work = work(j);
                factor(j, k) = cosine * old_factor + sine * old_work;
                work(j) = cosine * old_work - sine * old_factor;
            }
        } else {
            if (entry == Scalar(0)) {
                // The hyperbolic rotation is the identity, so a zero pivot
                // of a semidefinite factor can be passed through unchanged.
                continue;
            }
            if (!(std::abs(diagonal) > std::abs(entry))) {
                return k;
            }
            const Scalar sine = entry / diagonal;
            const Scalar magnitude = std::sqrt((Scalar(1) - sine) * (Scalar(1) + sine));
            if (!(magnitude > Scalar(0))) {
                return k;
            }
            const Scalar cosine = std::copysign(magnitude, diagonal);
            factor(k, k) = diagonal * cosine;
            for (Index j = k + 1; j < factor.rows(); ++j) {
                factor(j, k) = (factor(j, k) - sine * work(j)) / cosine;
                work(j) = cosine * work(j) - sine * factor(j, k);
            }
        }
    }
    return Index(-1);
}

/* For an upper factor U, update U^T U through its lower transpose view. */
template <int UpLo, typename MatrixType, typename VectorType>
typename MatrixType::Index rank_update(MatrixType& factor,
        const VectorType& vector, typename MatrixType::Scalar weight) {
    static_assert(UpLo == Eigen::Lower || UpLo == Eigen::Upper,
        "UKF rank updates require a lower or upper triangular factor");
    if constexpr (UpLo == Eigen::Lower) {
        return rank_update_lower(factor, vector, weight);
    } else {
        auto lower = factor.transpose();
        return rank_update_lower(lower, vector, weight);
    }
}

}
}
#endif
