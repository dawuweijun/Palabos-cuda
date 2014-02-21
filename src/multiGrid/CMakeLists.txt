set(MULTIGRID_SRC
	domainDivision3D.cpp
	multiGridDataProcessorWrapper2D.cpp
	multiGridManagement3D.cpp
	multiScale.cpp
	multiGrid2D.cpp
	multiGridDataProcessorWrapper3D.cpp
	multiGridOperations2D.cpp
	parallelizer2D.cpp
	multiGrid3D.cpp
	multiGridManagement2D.cpp
	multiGridOperations3D.cpp
	parallelizer3D.cpp
)
add_library(multiGrid STATIC ${MULTIGRID_SRC})
