/*
GCC 13 and newer can report array bounds in Eigen 3.2.8/3.4.0's optimized
fixed/bounded scalar QR, LLT and packet paths. They also report possibly
uninitialized values inside Eigen: matrixV in Eigen 3.4.0's quaternion
FromTwoVectors SVD (GCC 13), even though ComputeFullV initializes it for valid
input, and the eigenvectors of a bounded 1x1 SelfAdjointEigenSolver (GCC 15). Apply the diagnostic state while
parsing Eigen, then restore it before any UKF code. This header is force-included only by
project targets for that compiler/dependency combination; it is not exported
through the public CMake targets.
*/
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Warray-bounds"
#if __GNUC__ >= 13
#pragma GCC diagnostic ignored "-Wmaybe-uninitialized"
#endif
#include <Eigen/Dense>
#pragma GCC diagnostic pop
