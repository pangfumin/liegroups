CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

PROJECT(liegroup)
SET(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)

set(CMAKE_BUILD_TYPE DEBUG)

FIND_PACKAGE(OpenCV REQUIRED)
#set(OpenCV_INCLUDE_DIRS /home/pang/software/opencv3_install/include)
#set(OpenCV_LIBRARIES /home/pang/software/opencv3_install/lib)
include_directories(BEFORE ${OpenCV_INCLUDE_DIRS}) 

FIND_PACKAGE(Eigen3 REQUIRED)
include_directories(${EIGEN_INCLUDE_DIR}) 

FIND_PACKAGE(Ceres REQUIRED)
FIND_PACKAGE(Boost COMPONENTS program_options filesystem system thread chrono REQUIRED)
FIND_PACKAGE(OpenMP)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native -Wall -std=c++11 -fPIC")

SET(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3 -funroll-loops")
SET(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -Wall -Wno-unknown-pragmas -pedantic -g")

IF(OPENMP_FOUND)
    ADD_DEFINITIONS("-DOPENMP_FOUND")
    SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
ENDIF(OPENMP_FOUND)


include_directories(
	include) 


add_library(liegroup 
src/aff2.cpp
src/aff2_io.cpp
src/exp_coefs.cpp
src/matrix.cpp
src/se2.cpp
src/se2_io.cpp
src/se3.cpp
src/se3_io.cpp
src/sim2.cpp
src/sim2_io.cpp
src/sim3.cpp
src/sim3_io.cpp
src/sl3.cpp
src/sl3_io.cpp
src/so3.cpp
src/so3_io.cpp
src/spline.cpp
src/spline_io.cpp

)
#target_link_libraries(vins_app  VINS        )

add_subdirectory(test)





