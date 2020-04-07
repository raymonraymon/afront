PROJECT(Viewer)

IF (NOT DEFINED CMAKE_CXX_FLAGS)
   SET (CMAKE_CXX_FLAGS "-O2 -DREAL_IS_FLOAT -DNO_EXCEPTIONS -D_REENTRANT -fno-exceptions")
ENDIF (NOT DEFINED CMAKE_CXX_FLAGS)

FIND_PATH(GTB_INCLUDE_DIR gtb/gtb.hpp ../ ../.. /usr/include /usr/local/include ${HOME}/local/include)
INCLUDE_DIRECTORIES(${GTB_INCLUDE_DIR} ../)
ADD_LIBRARY(viewer STATIC debug_draw.cpp trackball.cpp viewer.cpp)
