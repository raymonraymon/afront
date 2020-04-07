

SET (CMAKE_CXX_FLAGS "-DUSENR -DREAL_IS_FLOAT -DNO_EXCEPTIONS -D_REENTRANT -fno-exceptions")

INCLUDE_DIRECTORIES(.. ../lib)


ADD_SUBDIRECTORY(../gtb)
ADD_SUBDIRECTORY(../lib/viewer)
ADD_SUBDIRECTORY(../lib/mlslib)
ADD_SUBDIRECTORY(../lib/ThreadLib)



PROJECT(afront)
LINK_DIRECTORIES(../gtb ../lib/mlslib ../lib/viewer ../lib/ThreadLib)



ADD_EXECUTABLE(afront
	crease.cpp      main.cpp                   output_controller_reeb.cpp  triangulate_mesh.cpp
	front.cpp       MarchingCubes.cpp          rg.cpp
	generaldef.cpp  output_controller_gui.cpp  triangulate_csg.cpp         triangulate_mls_smooth.cpp
	guidance.cpp    output_controller_hhm.cpp  triangulate_tet.cpp
	lsqr.cpp        output_controller_obj.cpp  triangulate_iso.cpp         triangulator.cpp
	edgeflipper.cpp FLF_io.cpp                 PC_io.cpp)


TARGET_LINK_LIBRARIES(afront viewer ThreadLib gtb mlslib
			     GL GLU glut m pthread
			     z)

