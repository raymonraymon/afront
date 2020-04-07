PROJECT(mlslib)

IF (NOT DEFINED CMAKE_CXX_FLAGS)
   SET(CMAKE_CXX_FLAGS "-O2 -DREAL_IS_FLOAT -DNO_EXCEPTIONS -D_REENTRANT -fno-exceptions")
ENDIF (NOT DEFINED CMAKE_CXX_FLAGS)

FIND_PATH(GTB_INCLUDE_DIR gtb/gtb.hpp ../.. /usr/include /usr/local/include ${HOME}/local/include)
FIND_PATH(VGC_STAT_INCLUDE_DIR stat/stat.h .. ../lib /usr/include /usr/local/include ${HOME}/local/include)

INCLUDE_DIRECTORIES(${GTB_INCLUDE_DIR} ${VGC_STAT_INCLUDE_DIR})

ADD_LIBRARY(mlslib STATIC CProjection.cpp knn.cpp NR/nrutil.cpp poly1.cpp poly2.cpp)
