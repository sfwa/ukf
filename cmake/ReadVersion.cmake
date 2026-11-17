# Reads the library version from its public version header, which is the one
# place it is written, so that the CMake package version cannot drift from it.
#
#   read_version(<header> <PREFIX> <version-var> <prerelease-var>)
#
# reads PREFIX_VERSION_MAJOR, _MINOR, _PATCH, _PRERELEASE and _STRING, sets
# <version-var> to MAJOR.MINOR.PATCH and <prerelease-var> to the pre-release
# (empty for a release), and fails the configure if _STRING is not
# MAJOR.MINOR.PATCH followed by -PRERELEASE when there is one.
function(read_version header prefix version_var prerelease_var)
    file(STRINGS "${header}" lines REGEX "^#define ${prefix}_VERSION_")
    foreach(part MAJOR MINOR PATCH)
        if(NOT lines MATCHES "#define ${prefix}_VERSION_${part} ([0-9]+)")
            message(FATAL_ERROR "${header} does not define ${prefix}_VERSION_${part}")
        endif()
        set(${part} ${CMAKE_MATCH_1})
    endforeach()
    foreach(part PRERELEASE STRING)
        if(NOT lines MATCHES "#define ${prefix}_VERSION_${part} \"([0-9A-Za-z.+-]*)\"")
            message(FATAL_ERROR "${header} does not define ${prefix}_VERSION_${part} as a string")
        endif()
        set(${part} "${CMAKE_MATCH_1}")
    endforeach()
    set(version ${MAJOR}.${MINOR}.${PATCH})
    set(expected ${version})
    if(NOT PRERELEASE STREQUAL "")
        if(NOT PRERELEASE MATCHES "^[0-9A-Za-z-]+(\\.[0-9A-Za-z-]+)*$")
            message(FATAL_ERROR "${prefix}_VERSION_PRERELEASE \"${PRERELEASE}\" is not a SemVer pre-release")
        endif()
        set(expected ${version}-${PRERELEASE})
    endif()
    if(NOT STRING STREQUAL expected)
        message(FATAL_ERROR "${header}: ${prefix}_VERSION_STRING is \"${STRING}\" but the numbers make \"${expected}\"")
    endif()
    set(${version_var} ${version} PARENT_SCOPE)
    set(${prerelease_var} "${PRERELEASE}" PARENT_SCOPE)
endfunction()
