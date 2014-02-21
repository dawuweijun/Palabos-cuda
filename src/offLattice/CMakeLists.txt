set(OFFLATTICE_SRC
      domainClustering3D.cpp
      marchingCube.cpp
)
add_library(offLattice STATIC ${OFFLATTICE_SRC})
