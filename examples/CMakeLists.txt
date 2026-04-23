cmake_minimum_required(VERSION 3.16)

# Solar System Example
add_executable(solar_system solar_system/main.cpp)
target_link_libraries(solar_system lagrange)
target_include_directories(solar_system PRIVATE ${CMAKE_SOURCE_DIR}/include)

# Rigid Body Example  
add_executable(rigid_body rigid_body/main.cpp)
target_link_libraries(rigid_body lagrange)
target_include_directories(rigid_body PRIVATE ${CMAKE_SOURCE_DIR}/include)

# Particle System Example
add_executable(particle_system particle_system/main.cpp)
target_link_libraries(particle_system lagrange)
target_include_directories(particle_system PRIVATE ${CMAKE_SOURCE_DIR}/include)
