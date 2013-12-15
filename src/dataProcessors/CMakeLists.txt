set(DATAPROCESSORS_SRC
metaStuffFunctional2D.cpp
metaStuffFunctional3D.cpp
metaStuffWrapper2D.cpp
metaStuffWrapper3D.cpp
)

add_library(dataProcessors SHARED ${DATAPROCESSORS_SRC})
