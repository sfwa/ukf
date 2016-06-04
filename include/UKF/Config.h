/* Selects the floating-point precision to use. */
#if defined(UKF_DOUBLE_PRECISION)
  typedef double real_t;
#elif defined(UKF_SINGLE_PRECISION)
  typedef float real_t;
#else
  #error "Define the floating-point precision using either UKF_DOUBLE_PRECISION or UKF_SINGLE_PRECISION"
#endif