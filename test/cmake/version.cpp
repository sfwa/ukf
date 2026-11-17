// UKF/Version.h stands alone: it needs no precision macro and no Eigen.
#include "UKF/Version.h"

#ifdef UKF_SINGLE_PRECISION
#error "This consumer must not choose a precision"
#endif
static_assert(UKF_VERSION_MAJOR >= 1, "UKF_VERSION_MAJOR must be defined");
static_assert(sizeof(UKF_VERSION_STRING) > 5, "UKF_VERSION_STRING must be defined");

int main() { return UKF_VERSION > 0 ? 0 : 1; }
