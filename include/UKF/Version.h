/*
Part of UKF. Distributed under the MIT licence: see LICENSE in the root of
the repository.
*/

#ifndef UKF_VERSION_H
#define UKF_VERSION_H

/*
The library version, under Semantic Versioning (https://semver.org).
RELEASING.md says what each part of the version promises. This header needs
no precision macro, so it may be included on its own, and UKF/Config.h (which
the filter headers include) includes it.

This header is the one place the version is written. CMake reads it for the
ukf package version (find_package(ukf 1.0)), and a configure fails if the
string below disagrees with the numbers. UKF_VERSION_PRERELEASE is empty in a
release and otherwise names the pre-release, such as "rc.1". Between releases
master carries the next planned version with the pre-release "dev".
*/
#define UKF_VERSION_MAJOR 1
#define UKF_VERSION_MINOR 0
#define UKF_VERSION_PATCH 0
#define UKF_VERSION_PRERELEASE "rc.1"
#define UKF_VERSION_STRING "1.0.0-rc.1"

/* One number for comparisons in the preprocessor: 10000 for 1.0.0, 10203 for
   1.2.3. */
#define UKF_VERSION \
  (UKF_VERSION_MAJOR * 10000 + UKF_VERSION_MINOR * 100 + UKF_VERSION_PATCH)

#endif
