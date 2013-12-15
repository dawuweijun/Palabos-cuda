set(IO_SRC
colormaps.cpp
mpiParallelIO.cpp
multiBlockReader2D.cpp
multiBlockReader3D.cpp
multiBlockWriter2D.cpp
multiBlockWriter3D.cpp
parallelIO.cpp
plbFiles.cpp
serializerIO.cpp
serializerIO_2D.cpp
serializerIO_3D.cpp
vtkDataOutput.cpp
vtkStructuredDataOutput.cpp
)

add_library(io SHARED ${IO_SRC})
