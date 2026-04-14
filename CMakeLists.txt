cmake_minimum_required(VERSION 3.18)
project(kokkos-mpfun
    VERSION 0.1.0
    LANGUAGES CXX
    DESCRIPTION "Kokkos port of MPFUN2020 arbitrary-precision floating-point library"
)

# ============================================================================
# Options
# ============================================================================
option(MPFUN_BUILD_TESTS "Build unit tests" ON)
option(MPFUN_BUILD_BENCHMARKS "Build benchmarks" ON)
option(MPFUN_BUILD_EXAMPLES "Build example programs" ON)
option(MPFUN_ENABLE_ASSERTIONS "Enable runtime assertions" ON)

# ============================================================================
# C++ Standard
# ============================================================================
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# ============================================================================
# Find Kokkos
# ============================================================================
find_package(Kokkos REQUIRED)
message(STATUS "Found Kokkos: ${Kokkos_DIR}")
message(STATUS "Kokkos version: ${Kokkos_VERSION}")
message(STATUS "Kokkos devices: ${Kokkos_DEVICES}")

# ============================================================================
# Header-only library target
# ============================================================================
add_library(mpfun INTERFACE)
add_library(mpfun::mpfun ALIAS mpfun)

target_include_directories(mpfun INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
)

target_link_libraries(mpfun INTERFACE Kokkos::kokkos)

target_compile_features(mpfun INTERFACE cxx_std_17)

if(MPFUN_ENABLE_ASSERTIONS)
    target_compile_definitions(mpfun INTERFACE MPFUN_ENABLE_ASSERTIONS=1)
endif()

# ============================================================================
# Tests
# ============================================================================
if(MPFUN_BUILD_TESTS)
    enable_testing()
    add_subdirectory(tests)
endif()

# ============================================================================
# Benchmarks
# ============================================================================
if(MPFUN_BUILD_BENCHMARKS)
    add_subdirectory(benchmarks)
endif()

# ============================================================================
# Examples
# ============================================================================
if(MPFUN_BUILD_EXAMPLES)
    add_subdirectory(examples)
endif()

# ============================================================================
# Install
# ============================================================================
include(GNUInstallDirs)
include(CMakePackageConfigHelpers)

install(TARGETS mpfun
    EXPORT mpfun-targets
    INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)

install(DIRECTORY include/
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)

install(EXPORT mpfun-targets
    FILE mpfun-targets.cmake
    NAMESPACE mpfun::
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/mpfun
)

configure_package_config_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/cmake/mpfun-config.cmake.in
    ${CMAKE_CURRENT_BINARY_DIR}/mpfun-config.cmake
    INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/mpfun
)

write_basic_package_version_file(
    ${CMAKE_CURRENT_BINARY_DIR}/mpfun-config-version.cmake
    VERSION ${PROJECT_VERSION}
    COMPATIBILITY SameMajorVersion
)

install(FILES
    ${CMAKE_CURRENT_BINARY_DIR}/mpfun-config.cmake
    ${CMAKE_CURRENT_BINARY_DIR}/mpfun-config-version.cmake
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/mpfun
)
