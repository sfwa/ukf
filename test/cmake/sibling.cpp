#if defined(UKF_SINGLE_PRECISION) || defined(UKF_DOUBLE_PRECISION)
#error "UKF precision leaked to an unrelated target"
#endif
#ifdef _MSVC_LANG
static_assert(_MSVC_LANG == 201402L, "UKF changed the parent's C++ standard");
#else
static_assert(__cplusplus == 201402L, "UKF changed the parent's C++ standard");
#endif

// This should compile even when UKF's own targets enable -Wall -Wextra -Werror.
int main(int argc, char**) { return 0; }
