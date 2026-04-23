cmake_minimum_required(VERSION 3.20)
project(Lagrange VERSION 1.0.0 LANGUAGES C)

# ─────────────────────────────────────────────────────────────────────────────
# C11 – keep it modern
# ─────────────────────────────────────────────────────────────────────────────
set(CMAKE_C_STANDARD 11)
set(CMAKE_C_STANDARD_REQUIRED ON)
set(CMAKE_C_EXTENSIONS OFF)

# Default to Release so naive builds get reasonable performance
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build." FORCE)
    set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif()

# ─────────────────────────────────────────────────────────────────────────────
# Options
# ─────────────────────────────────────────────────────────────────────────────
option(LAGRANGE_BUILD_EXAMPLES          "Build example applications" ON)
option(LAGRANGE_BUILD_TESTS             "Build tests" ON)
option(LAGRANGE_ENABLE_NATIVE_OPTIMIZATIONS "Enable host-native CPU optimizations" ON)
option(LAGRANGE_ENABLE_LTO              "Enable Link-Time Optimization" ON)
option(LAGRANGE_ENABLE_SANITIZERS       "Enable ASan/UBSan on Debug builds" OFF)
option(LAGRANGE_WITH_CSPICE             "Fetch NASA CSPICE for astro validation harness" OFF)

# ─────────────────────────────────────────────────────────────────────────────
# Public interface library – only what consumers need
# ─────────────────────────────────────────────────────────────────────────────
add_library(lagrange INTERFACE)
target_include_directories(lagrange INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
)
target_compile_features(lagrange INTERFACE c_std_11)

# ─────────────────────────────────────────────────────────────────────────────
# Internal build options – warnings, sanitizers, native tuning.
# Kept separate so they do NOT leak to users via find_package(lagrange).
# ─────────────────────────────────────────────────────────────────────────────
add_library(lagrange_internal INTERFACE)
target_link_libraries(lagrange_internal INTERFACE lagrange)

if(LAGRANGE_ENABLE_NATIVE_OPTIMIZATIONS)
    include(CheckCCompilerFlag)
    if(CMAKE_C_COMPILER_ID MATCHES "GNU|Clang|AppleClang")
        check_c_compiler_flag("-march=native" COMPILER_SUPPORTS_MARCH_NATIVE)
        if(COMPILER_SUPPORTS_MARCH_NATIVE)
            target_compile_options(lagrange_internal INTERFACE -march=native -mtune=native)
        elseif(CMAKE_SYSTEM_PROCESSOR MATCHES "aarch64|arm64|ARM64")
            target_compile_options(lagrange_internal INTERFACE -mcpu=native -mtune=native)
        elseif(CMAKE_SYSTEM_PROCESSOR MATCHES "x86_64|AMD64")
            target_compile_options(lagrange_internal INTERFACE -march=x86-64-v3)
        endif()
    elseif(MSVC)
        target_compile_options(lagrange_internal INTERFACE /arch:AVX2)
    endif()
endif()

if(CMAKE_C_COMPILER_ID MATCHES "GNU|Clang|AppleClang")
    target_compile_options(lagrange_internal INTERFACE
        -Wall
        -Wextra
        -Wpedantic
        -Werror
        -Wshadow
        -Wstrict-prototypes
        -Wmissing-prototypes
        -Wstrict-aliasing=2
        -Wfloat-equal
        -Wundef
        -Wcast-align
        -Wno-unused-function
        -Wno-unused-parameter
        -Wno-missing-prototypes
        # -Wconversion disabled temporarily: legacy code has extensive sign/int-float
        # conversions that need careful cleanup. Re-enable after fixing koopman.h,
        # fractional.h, and wavelet.h.
    )
elseif(MSVC)
    target_compile_options(lagrange_internal INTERFACE /W4 /WX /permissive-)
endif()

if(LAGRANGE_ENABLE_SANITIZERS AND CMAKE_BUILD_TYPE STREQUAL "Debug")
    if(CMAKE_C_COMPILER_ID MATCHES "GNU|Clang|AppleClang")
        target_compile_options(lagrange_internal INTERFACE -fsanitize=address,undefined -fno-omit-frame-pointer)
        target_link_options(lagrange_internal INTERFACE -fsanitize=address,undefined)
    endif()
endif()

if(LAGRANGE_ENABLE_LTO)
    include(CheckIPOSupported)
    check_ipo_supported(RESULT ipo_supported OUTPUT ipo_error)
    if(ipo_supported)
        set(CMAKE_INTERPROCEDURAL_OPTIMIZATION ON)
    else()
        message(WARNING "IPO/LTO not supported by this compiler: ${ipo_error}")
    endif()
endif()

# ─────────────────────────────────────────────────────────────────────────────
# FetchContent – greatest.h test harness (git clone, tiny repo)
# ─────────────────────────────────────────────────────────────────────────────
include(FetchContent)

FetchContent_Declare(
    greatest
    GIT_REPOSITORY https://github.com/silentbicycle/greatest.git
    GIT_TAG        release
    GIT_SHALLOW    TRUE
)
FetchContent_MakeAvailable(greatest)

# ─────────────────────────────────────────────────────────────────────────────
# CSPICE integration – NASA gold-standard ephemeris oracle
# Only pulled when LAGRANGE_WITH_CSPICE=ON and tests are enabled.
# ─────────────────────────────────────────────────────────────────────────────
if(LAGRANGE_WITH_CSPICE AND LAGRANGE_BUILD_TESTS)
    set(CSPICE_URL "")
    if(APPLE)
        if(CMAKE_SYSTEM_PROCESSOR MATCHES "arm64|aarch64|ARM64")
            set(CSPICE_URL "https://naif.jpl.nasa.gov/pub/naif/toolkit/C/MacM1_OSX_clang_64bit/packages/cspice.tar.Z")
        else()
            set(CSPICE_URL "https://naif.jpl.nasa.gov/pub/naif/toolkit/C/MacIntel_OSX_AppleC_64bit/packages/cspice.tar.Z")
        endif()
    elseif(UNIX)
        if(CMAKE_SYSTEM_PROCESSOR MATCHES "x86_64|AMD64")
            set(CSPICE_URL "https://naif.jpl.nasa.gov/pub/naif/toolkit/C/PC_Linux_GCC_64bit/packages/cspice.tar.Z")
        else()
            message(FATAL_ERROR "CSPICE prebuilt not available for ${CMAKE_SYSTEM_NAME}/${CMAKE_SYSTEM_PROCESSOR}. NAIF only provides x86_64 Linux prebuilds. Set LAGRANGE_WITH_CSPICE=OFF or build CSPICE from source.")
        endif()
    elseif(WIN32)
        set(CSPICE_URL "https://naif.jpl.nasa.gov/pub/naif/toolkit/C/PC_Windows_VisualC_64bit/packages/cspice.zip")
    else()
        message(FATAL_ERROR "CSPICE prebuilt not available for ${CMAKE_SYSTEM_NAME}/${CMAKE_SYSTEM_PROCESSOR}")
    endif()

    FetchContent_Declare(cspice URL ${CSPICE_URL})
    FetchContent_MakeAvailable(cspice)

    # CSPICE archives extract with a top-level "cspice/" directory.
    # Handle both nested and flat extraction gracefully.
    set(CSPICE_ROOT "${cspice_SOURCE_DIR}/cspice")
    if(NOT EXISTS "${CSPICE_ROOT}")
        set(CSPICE_ROOT "${cspice_SOURCE_DIR}")
    endif()

    if(WIN32)
        set(CSPICE_LIB "${CSPICE_ROOT}/lib/cspice.lib")
    else()
        set(CSPICE_LIB "${CSPICE_ROOT}/lib/cspice.a")
    endif()

    if(NOT EXISTS "${CSPICE_LIB}")
        message(FATAL_ERROR "CSPICE library not found at ${CSPICE_LIB}")
    endif()

    add_library(cspice STATIC IMPORTED GLOBAL)
    set_target_properties(cspice PROPERTIES
        IMPORTED_LOCATION "${CSPICE_LIB}"
        INTERFACE_INCLUDE_DIRECTORIES "${CSPICE_ROOT}/include"
    )
    if(UNIX OR APPLE)
        target_link_libraries(cspice INTERFACE m)
    endif()

    message(STATUS "CSPICE integrated: ${CSPICE_ROOT}")
endif()

# ─────────────────────────────────────────────────────────────────────────────
# Math library – only where it actually exists
# ─────────────────────────────────────────────────────────────────────────────
set(MATH_LIB "")
if(UNIX OR APPLE)
    set(MATH_LIB m)
endif()

# ─────────────────────────────────────────────────────────────────────────────
# Tests
# ─────────────────────────────────────────────────────────────────────────────
if(LAGRANGE_BUILD_TESTS)
    enable_testing()

    function(add_lagrange_test name source_file)
        add_executable(${name} ${source_file})
        target_link_libraries(${name} PRIVATE lagrange_internal ${MATH_LIB})
        target_include_directories(${name} PRIVATE ${greatest_SOURCE_DIR})
        add_test(NAME ${name} COMMAND ${name})
    endfunction()

    add_lagrange_test(test_math      tests/test_math.c)
    add_lagrange_test(test_physics   tests/test_physics.c)
    add_lagrange_test(test_koopman   tests/test_koopman.c)
    add_lagrange_test(test_wavelet   tests/test_wavelet.c)
    add_lagrange_test(test_orbital   tests/test_orbital.c)
    add_lagrange_test(test_attitude  tests/test_attitude.c)
    add_lagrange_test(test_gravity   tests/test_gravity.c)
    add_lagrange_test(test_rings     tests/test_rings.c)

    if(LAGRANGE_WITH_CSPICE AND TARGET cspice)
        if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/tests/test_astro_validation.c")
            add_executable(test_astro_validation tests/test_astro_validation.c)
            target_link_libraries(test_astro_validation PRIVATE lagrange_internal cspice ${MATH_LIB})
            target_include_directories(test_astro_validation PRIVATE ${greatest_SOURCE_DIR})
            add_test(NAME test_astro_validation COMMAND test_astro_validation)
            message(STATUS "CSPICE astro validation test enabled")
        else()
            message(WARNING "LAGRANGE_WITH_CSPICE=ON but tests/test_astro_validation.c not found. Create it to enable validation.")
        endif()
    endif()
endif()

# ─────────────────────────────────────────────────────────────────────────────
# Examples
# ─────────────────────────────────────────────────────────────────────────────
if(LAGRANGE_BUILD_EXAMPLES)
    function(add_lagrange_example name source_file)
        add_executable(${name} ${source_file})
        target_link_libraries(${name} PRIVATE lagrange_internal ${MATH_LIB})
    endfunction()

    add_lagrange_example(bouncing_ball examples/bouncing_ball.c)
    add_lagrange_example(orbit         examples/orbit.c)
    add_lagrange_example(stacking      examples/stacking.c)
endif()

# ─────────────────────────────────────────────────────────────────────────────
# Install / export – find_package(lagrange) ready
# ─────────────────────────────────────────────────────────────────────────────
include(GNUInstallDirs)

install(TARGETS lagrange
    EXPORT lagrange-targets
)

install(EXPORT lagrange-targets
    FILE lagrange-targets.cmake
    NAMESPACE lagrange::
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/lagrange
)

install(DIRECTORY include/
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)

include(CMakePackageConfigHelpers)
write_basic_package_version_file(
    "${CMAKE_CURRENT_BINARY_DIR}/lagrange-config-version.cmake"
    VERSION ${PROJECT_VERSION}
    COMPATIBILITY SameMajorVersion
)

configure_package_config_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/cmake/lagrange-config.cmake.in"
    "${CMAKE_CURRENT_BINARY_DIR}/lagrange-config.cmake"
    INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/lagrange
)

install(FILES
    "${CMAKE_CURRENT_BINARY_DIR}/lagrange-config.cmake"
    "${CMAKE_CURRENT_BINARY_DIR}/lagrange-config-version.cmake"
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/lagrange
)
