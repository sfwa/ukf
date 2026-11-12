// No Eigen include: UKF headers must supply their own geometry dependencies.
#include "UKF/Core.h"
#include <type_traits>

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

// Deliberately unused: the library must not export its warning policy.
int main(int argc, char**) {
    const UKF::Quaternion orientation = UKF::Quaternion::Identity();
    const UKF::Vector<3> vector = UKF::Vector<3>::Ones();
    return (orientation * vector - vector).squaredNorm() == real_t(0) ? 0 : 1;
}
