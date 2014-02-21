set(PARALLELISM_SRC
	mpiManager.cpp
	parallelBlockCommunicator2D.cpp
	parallelBlockCommunicator3D.cpp
	parallelStatistics.cpp
	sendRecvPool.cpp
)

add_library(parallelism STATIC ${PARALLELISM_SRC})
