file(GLOB_RECURSE CPP_FILE_LIST *.cpp)

# message(FATAL_ERROR "{${CPP_FILE_LIST}")

add_definitions(-DPLB_MPI_PARALLEL -DPLB_USE_POSIX -DPLB_DEBUG)

add_library(Palabos SHARED ${CPP_FILE_LIST})

target_link_libraries(Palabos tinyxml ${MPI_CXX_LIBRARIES})

install(TARGETS Palabos DESTINATION ${INSTALL_DIR_LIB} COMPONENT library)
