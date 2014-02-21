cmake_minimum_required(VERSION 2.6)
project(Palabos-Cuda)

include_directories(/usr/lib/openmpi/include/)


#######################################################################
#
#######################################################################
find_package(Eigen3 REQUIRED)
if(EIGEN3_FOUND)
	include_directories(${EIGEN3_INCLUDE_DIR})
elseif(EIGEN3_FOUND)
	message(FATAL_ERROR "Eigen3 Not Found!")
endif(EIGEN3_FOUND)

#######################################################################
#
#######################################################################
find_package(PkgConfig)
if(PKG_CONFIG_FOUND)
      pkg_check_modules(TINYXML REQUIRED tinyxml)
      include_directories(TINYXML_INCLUDE_DIRS)
elseif(PKG_CONFIG_FOUND)
      message(FATAL_ERROR "Pkg-config Not Found!")
endif(PKG_CONFIG_FOUND)

#######################################################################
#
#######################################################################
# find_package(CUDA)
# if(CUDA_FOUND)
# 	include_directories(${CUDA_INCLUDE_DIRS})
# 	message(FATAL_ERROR ${CUDA_INCLUDE_DIRS})
# elseif(CUDA_FOUND)
# 	message(FATAL_ERROR "Cuda Not Found!")
# endif(CUDA_FOUND)

include_directories(/usr/local/cuda-6.0/include/)

#######################################################################
#
#######################################################################
find_package(MPI)
if(MPI_CXX_FOUND)
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${MPI_CXX_COMPILE_FLAGS}")
# 	set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS}:${OpenMP_C_FLAGS}")
	include_directories(${MPI_CXX_INCLUDE_PATH})
	add_definitions(-DPLB_MPI_PARALLEL)
elseif(MPI_CXX_FOUND)
	message("MPI Not Support")
endif(MPI_CXX_FOUND)

#######################################################################
#
#######################################################################
if("${CMAKE_C_COMPILER_ID}" STREQUAL "GNU")
    add_definitions(-DPLB_USE_GCC)
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")

elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "BGP")
    add_definitions(-DPLB_USE_BGP)
endif()

include_directories(include)

add_subdirectory(src)
add_subdirectory(include)
add_subdirectory(utility)
add_subdirectory(examples)



#######################################################################
# uninstall target
#######################################################################

configure_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/cmake_uninstall.cmake.in"
    "${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake"
    IMMEDIATE @ONLY)

add_custom_target(uninstall
    COMMAND ${CMAKE_COMMAND} -P ${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake)
