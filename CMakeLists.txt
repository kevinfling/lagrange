cmake_minimum_required(VERSION 3.16)
project(Lagrange VERSION 1.0.0 LANGUAGES C)

set(CMAKE_C_STANDARD 11)
set(CMAKE_C_STANDARD_REQUIRED ON)

# Options
option(LAGRANGE_BUILD_EXAMPLES "Build example applications" ON)
option(LAGRANGE_BUILD_TESTS "Build tests" ON)

# Header-only interface library
add_library(lagrange INTERFACE)
target_include_directories(lagrange INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
)

# Installation
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

# Tests
if(LAGRANGE_BUILD_TESTS)
    enable_testing()
    
    add_executable(test_math tests/test_math.c)
    target_link_libraries(test_math lagrange m)
    add_test(NAME math COMMAND test_math)
    
    add_executable(test_physics tests/test_physics.c)
    target_link_libraries(test_physics lagrange m)
    add_test(NAME physics COMMAND test_physics)
    
    add_executable(test_koopman tests/test_koopman.c)
    target_link_libraries(test_koopman lagrange m)
    add_test(NAME koopman COMMAND test_koopman)
    
    add_executable(test_wavelet tests/test_wavelet.c)
    target_link_libraries(test_wavelet lagrange m)
    add_test(NAME wavelet COMMAND test_wavelet)
    
    add_executable(test_orbital tests/test_orbital.c)
    target_link_libraries(test_orbital lagrange m)
    add_test(NAME orbital COMMAND test_orbital)
    
    add_executable(test_attitude tests/test_attitude.c)
    target_link_libraries(test_attitude lagrange m)
    add_test(NAME attitude COMMAND test_attitude)
    
    add_executable(test_gravity tests/test_gravity.c)
    target_link_libraries(test_gravity lagrange m)
    add_test(NAME gravity COMMAND test_gravity)
    
    add_executable(test_rings tests/test_rings.c)
    target_link_libraries(test_rings lagrange m)
    add_test(NAME rings COMMAND test_rings)
endif()

# Examples
if(LAGRANGE_BUILD_EXAMPLES)
    add_executable(bouncing_ball examples/bouncing_ball.c)
    target_link_libraries(bouncing_ball lagrange m)
    
    add_executable(orbit examples/orbit.c)
    target_link_libraries(orbit lagrange m)
    
    add_executable(stacking examples/stacking.c)
    target_link_libraries(stacking lagrange m)
endif()
