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
find_package(CUDA)
if(CUDA_FOUND)
	include_directories(${CUDA_INCLUDE_DIRS})
elseif(CUDA_FOUND)
	message(FATAL_ERROR "Cuda Not Found!")
endif(CUDA_FOUND)

add_subdirectory(src)
add_subdirectory(include)
add_subdirectory(utility)
add_subdirectory(examples)