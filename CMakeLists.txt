cmake_minimum_required (VERSION 3.0)

project(MultiSensorAutoCalib CXX)

if (MSVC)
    set(CMAKE_CXX_FLAGS "-std=c++11")
    set(src src)
    set(EIGEN_INCLUDE_DIR D:/eigen3)

    find_package(Eigen3 REQUIRED)

    include_directories(${EIGEN_INCLUDE_DIR})
    include_directories(${src})

    add_executable(${PROJECT_NAME}
                   ${src}/main.cxx
                   )
else()

endif()
