cmake_minimum_required(VERSION 3.24)
project(embedded_methods_cpp)
find_package(Eigen3 3.4 REQUIRED NO_MODULE)
set(CMAKE_CXX_STANDARD 17)

add_executable(embedded_methods_cpp src/main.cpp src/EmbeddedRungeKutta.h src/RK45Solver.h src/RKUtils.h src/SRK64Solver.h src/RKLogger.h src/RK54Solver.h src/DOPRI54Solver.h)
target_link_libraries(embedded_methods_cpp Eigen3::Eigen)
