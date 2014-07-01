add_executable(tutorial_2_1 tutorial_2_1.cpp)

target_link_libraries(tutorial_2_1 ${MPI_CXX_LIBRARIES} Palabos tinyxml)

# install(TARGETS tutorial_1_1 DESTINATION ${CMAKE_INSTALL_PREFIX}/examples)
