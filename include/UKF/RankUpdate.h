/* SPDX-License-Identifier: MIT */
#ifndef UKF_RANK_UPDATE_H
#define UKF_RANK_UPDATE_H

#include <Eigen/Core>
#include <cmath>
#include <limits>
#include <type_traits>

namespace UKF {
namespace Detail {

/*
Update L L^T by weight * V V^T, where each column of V is one vector, using
orthogonal (positive weight) or hyperbolic (negative weight) rotations. Each
column of L is rotated against every vector in turn. These are the rotations
of one rank-one update after another, in an order that lets a processor
overlap them. Keep the workspace bounded, unlike Eigen's LLT rank-update
temporary in unpatched Eigen releases. Only the lower triangle is read or
changed. Return -1 on success, or the index of the first pivot a downdate
would make nonpositive; callers must discard a partially completed failed
downdate.

A pivot passes through a downdate unchanged if its update entry is zero, or
if the downdate would fail there but both are negligible beside the rest of
the pivot's row. Such a direction, as in a semidefinite factor, has no
variance of its own to remove beyond roundoff.
*/
template <typename MatrixType, typename VectorsType>
typename MatrixType::Index rank_update_lower(MatrixType& factor,
        const VectorsType& vectors, typename MatrixType::Scalar weight) {
    using Scalar = typename MatrixType::Scalar;
    using Index = typename MatrixType::Index;
    static_assert(std::is_floating_point<Scalar>::value,
        "UKF rank updates require real floating-point factors");
    static_assert(MatrixType::MaxRowsAtCompileTime != Eigen::Dynamic &&
        VectorsType::MaxColsAtCompileTime != Eigen::Dynamic,
        "UKF rank updates require fixed or bounded matrix storage");
    eigen_assert(factor.rows() == factor.cols() && factor.rows() == vectors.rows());
    if (weight == Scalar(0)) {
        return Index(-1);
    }

    const Scalar tolerance = std::sqrt(std::numeric_limits<Scalar>::epsilon());
    Eigen::Matrix<Scalar, MatrixType::RowsAtCompileTime, VectorsType::ColsAtCompileTime,
        MatrixType::MaxRowsAtCompileTime == 1 ? Eigen::RowMajor : Eigen::ColMajor,
        MatrixType::MaxRowsAtCompileTime, VectorsType::MaxColsAtCompileTime> work =
        std::sqrt(std::abs(weight)) * vectors;

    for (Index k = 0; k < factor.rows(); ++k) {
        for (Index i = 0; i < work.cols(); ++i) {
            const Scalar diagonal = factor(k, k);
            const Scalar entry = work(k, i);
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
                    const Scalar old_work = work(j, i);
                    factor(j, k) = cosine * old_factor + sine * old_work;
                    work(j, i) = cosine * old_work - sine * old_factor;
                }
            } else {
                const bool reducible = std::abs(diagonal) > std::abs(entry);
                if (entry == Scalar(0) || (!reducible &&
                        std::abs(entry) <= tolerance * factor.row(k).head(k).norm())) {
                    // The hyperbolic rotation is the identity, so a zero pivot
                    // of a semidefinite factor can be passed through. Still
                    // return a canonical nonnegative pivot, as below.
                    if (diagonal < Scalar(0)) {
                        for (Index j = k; j < factor.rows(); ++j) {
                            factor(j, k) = -factor(j, k);
                        }
                    }
                    continue;
                }
                if (!reducible) {
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
                    factor(j, k) = (factor(j, k) - sine * work(j, i)) / cosine;
                    work(j, i) = cosine * work(j, i) - sine * factor(j, k);
                }
            }
        }
    }
    return Index(-1);
}

/* For an upper factor U, update U^T U through its lower transpose view. */
template <int UpLo, typename MatrixType, typename VectorsType>
typename MatrixType::Index rank_update(MatrixType& factor,
        const VectorsType& vectors, typename MatrixType::Scalar weight) {
    static_assert(UpLo == Eigen::Lower || UpLo == Eigen::Upper,
        "UKF rank updates require a lower or upper triangular factor");
    if constexpr (UpLo == Eigen::Lower) {
        return rank_update_lower(factor, vectors, weight);
    } else {
        auto lower = factor.transpose();
        return rank_update_lower(lower, vectors, weight);
    }
}

}
}
#endif
