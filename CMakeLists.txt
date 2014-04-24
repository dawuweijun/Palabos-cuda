project(Palabos)
cmake_minimum_required(VERSION 2.8)


find_package(Eigen3 REQUIRED)

if(${EIGEN3_FOUND})
		include_directories(${EIGEN3_INCLUDE_DIR})
else(${EIGEN3_FOUND})
		message(FATAL_ERROR "Required Package Eigen3 Not Found")
endif(${EIGEN3_FOUND})


find_package(PkgConfig REQUIRED)

if(${PKG_CONFIG_FOUND})
	pkg_check_modules(TinyXML REQUIRED tinyxml)
	if(TinyXML_FOUND)
		include_directories(${TinyXML_INCLUDE_DIRS})
#  		message(FATAL_ERROR ${TinyXML_INCLUDE_DIRS})
	else(TinyXML_FOUND)
		message(FATAL_ERROR "Required Package TinyXML Not Found")
	endif(TinyXML_FOUND)
else(${PKG_CONFIG_FOUND})
	message(FATAL_ERROR "Required Package PkgConfig Not Found")
endif(${PKG_CONFIG_FOUND})

find_package(MPI REQUIRED)

if(${MPI_FOUND})
	include_directories(${MPI_CXX_INCLUDE_PATH})
	add_definitions(-DPLB_MPI_PARALLEL)
else(${MPI_FOUND})
	message(FATAL_ERROR "Required Package MPI Not Found")
endif(${MPI_FOUND})

include_directories(src)

set(ALLOW_BUILD_HEADERS 0 CACHE TYPE BOOL)
set(ALLOW_BUILD_PACKAGE 0 CACHE TYPE BOOL)

# SET RELEVANT VARS FOR PACKAGE
if (ALLOW_BUILD_PACKAGE)
	if (APPLE)
	elseif (UNIX)
	    SET(PACKAGE_TYPE "DEB")    # DEB or RPM
	endif (APPLE)

	SET(PACKAGE_MAINTAINER "HuiJie Zhang <yitianbuji@gmail.com>")
	#SET(PACKAGE_DESCRIPTION_FILE "${CMAKE_CURRENT_SOURCE_DIR}/about.txt")
	SET(PACKAGE_LICENSE_FILE "${CMAKE_CURRENT_SOURCE_DIR}/LICENSE")

	SET(PACKAGE_VERSION_MAJOR "1")
	SET(PACKAGE_VERSION_MINOR "4")
	SET(PACKAGE_VERSION_PATCH "1")

	SET(PACKAGE_SYSTEM_NAME "amd64")
endif (ALLOW_BUILD_PACKAGE)


if (APPLE)
    set(INSTALL_DIR_BIN "/usr/local/bin")
    set(INSTALL_DIR_LIB "/usr/local/lib")
    set(INSTALL_DIR_PLUGINS "/usr/local/bin")
    set(INSTALL_DIR_HEADERS "/usr/local/include")
elseif (UNIX)
    set(INSTALL_DIR_BIN "/usr/bin")
    set(INSTALL_DIR_LIB "/usr/lib")
    SET(INSTALL_DIR_PLUGINS "/usr/bin")
    set(INSTALL_DIR_HEADERS "/usr/include")
endif (APPLE)


# CREATE PACKAGE
if (ALLOW_BUILD_PACKAGE)
	if (APPLE)
	    # ...
	elseif (UNIX)
	    # COMMON
	    SET(CPACK_DEBIAN_PACKAGE_MAINTAINER ${PACKAGE_MAINTAINER})

	    SET(CPACK_PACKAGE_DESCRIPTION_FILE ${PACKAGE_DESCRIPTION_FILE})
	    SET(CPACK_RESOURCE_FILE_LICENSE ${PACKAGE_LICENSE_FILE})

	    SET(CPACK_PACKAGE_VERSION_MAJOR ${PACKAGE_VERSION_MAJOR})
	    SET(CPACK_PACKAGE_VERSION_MINOR ${PACKAGE_VERSION_MINOR})
	    SET(CPACK_PACKAGE_VERSION_PATCH ${PACKAGE_VERSION_PATCH})

	    SET(CPACK_SYSTEM_NAME ${PACKAGE_SYSTEM_NAME})

	    # DEB SPECIFIC
	    if (PACKAGE_TYPE STREQUAL "DEB")
	        SET(CPACK_GENERATOR "DEB")
	        SET(CPACK_DEBIAN_PACKAGE_SHLIBDEPS ON) # autogenerate dependency information (package d-shlibs should be installed)
	        #依赖Eigen3，tinyxml2,mpi,imagemagick,vtk6-dev
	        SET(CPACK_DEBIAN_PACKAGE_DEPENDS "vtk6 (>= 6.0.0), libboost-system1.54.0 (>= 1.54.0)")
	        #推荐paraview,
	    endif (PACKAGE_TYPE STREQUAL "DEB")

	    # RPM SPECIFIC
	    if (PACKAGE_TYPE STREQUAL "RPM")
	        SET(CPACK_GENERATOR "RPM")

	    endif(PACKAGE_TYPE STREQUAL "RPM")
	    INCLUDE(CPack)
	endif (APPLE)
endif (ALLOW_BUILD_PACKAGE)

add_definitions(-DPLB_CUDA_DISABLED)

include_directories(include)

add_subdirectory(utility)

add_subdirectory(examples)

add_subdirectory(src)

install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/include/
		DESTINATION ${INSTALL_DIR_HEADERS}/Palabos-1.4
		)
