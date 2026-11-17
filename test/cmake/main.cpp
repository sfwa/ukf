// No Eigen include: UKF headers must supply their own geometry dependencies.
#include <type_traits>

#include "UKF/Core.h"

#if defined(UKF_SINGLE_PRECISION) && defined(UKF_DOUBLE_PRECISION)
#error "A precision target must not define both precisions"
#endif
#ifdef _MSVC_LANG
static_assert(_MSVC_LANG >= 201703L, "UKF must propagate C++17 to consumers");
#else
static_assert(__cplusplus >= 201703L, "UKF must propagate C++17 to consumers");
#endif
#ifdef EXPECT_DOUBLE
static_assert(std::is_same_v<real_t, double>, "Expected double precision");
#else
static_assert(std::is_same_v<real_t, float>, "Expected single precision");
#endif

#ifdef EXPECT_EIGEN_UNROLLING_LIMIT
static_assert(EIGEN_UNROLLING_LIMIT == EXPECT_EIGEN_UNROLLING_LIMIT,
              "UKF_EIGEN_UNROLLING_LIMIT must reach consumers");
#endif

/* The header's version must be the CMake package's. */
constexpr bool same_string(const char* a, const char* b) {
  return *a == *b && (*a == '\0' || same_string(a + 1, b + 1));
}
static_assert(UKF_VERSION == UKF_VERSION_MAJOR * 10000 +
                                 UKF_VERSION_MINOR * 100 + UKF_VERSION_PATCH,
              "UKF_VERSION must be made from its parts");
#ifdef EXPECT_VERSION_STRING
static_assert(
    same_string(UKF_VERSION_STRING, EXPECT_VERSION_STRING),
    "UKF_VERSION_STRING must match the version of the ukf CMake package");
#endif

// Deliberately unused: the library must not export its warning policy.
int main(int argc, char**) {
  const UKF::Quaternion orientation = UKF::Quaternion::Identity();
  const UKF::Vector<3> vector = UKF::Vector<3>::Ones();
  return (orientation * vector - vector).squaredNorm() == real_t(0) ? 0 : 1;
}
