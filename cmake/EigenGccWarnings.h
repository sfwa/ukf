/*
GCC 13 and newer can report array bounds in Eigen 3.4.0's optimized fixed/bounded
scalar QR and packet paths. Apply the diagnostic state while parsing Eigen,
then restore it before any UKF code. This header is force-included only by
project targets for that compiler/dependency combination; it is not exported
through the public CMake targets.
*/
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Warray-bounds"
#include <Eigen/Dense>
#pragma GCC diagnostic pop
