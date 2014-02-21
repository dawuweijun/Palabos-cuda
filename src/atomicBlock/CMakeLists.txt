set(ATOMICBLOCK_SRC
	atomicBlock2D.cpp
	atomicBlock3D.cpp
	atomicBlockOperations2D.cpp
	atomicBlockOperations3D.cpp
	atomicBlockSerializer2D.cpp
	atomicBlockSerializer3D.cpp
	atomicContainerBlock2D.cpp
	atomicContainerBlock3D.cpp
	dataProcessingFunctional2D.cpp
	dataProcessingFunctional3D.cpp
	dataProcessor2D.cpp
	dataProcessor3D.cpp
	dataProcessorWrapper2D.cpp
	dataProcessorWrapper3D.cpp
	reductiveDataProcessingFunctional2D.cpp
	reductiveDataProcessingFunctional3D.cpp
	reductiveDataProcessorWrapper2D.cpp
	reductiveDataProcessorWrapper3D.cpp
)

add_library(atomicBlock STATIC ${ATOMICBLOCK_SRC})
