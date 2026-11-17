/*
GCC 13 and newer can report array bounds in Eigen 3.2.8/3.4.0's optimized
fixed/bounded scalar QR, LLT and packet paths. GCC 13 also reports a possibly
uninitialized matrixV in Eigen 3.4.0's quaternion FromTwoVectors SVD, even though
ComputeFullV initializes it for valid input. Apply the diagnostic state while
parsing Eigen, then restore it before any UKF code. This header is force-included only by
project targets for that compiler/dependency combination; it is not exported
through the public CMake targets.
*/
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Warray-bounds"
#if __GNUC__ == 13
#pragma GCC diagnostic ignored "-Wmaybe-uninitialized"
#endif
#include <Eigen/Dense>
#pragma GCC diagnostic pop
