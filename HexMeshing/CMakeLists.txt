# CMakeList.txt : CMake project for HexMeshing, include source and define
# project specific logic here.
#
cmake_minimum_required (VERSION 2.8)

find_package(VTK REQUIRED)
include(${VTK_USE_FILE})
include_directories(
	${VTK_INCLUDE_DIRS}
)
# Add source to this project's executable.
set(SRC_HEX_MESHING
	mesh.cpp
	mesh.h
)
add_executable (HexMeshing ${SRC_HEX_MESHING} "HexMeshing.cpp" "HexMeshing.h")

target_link_libraries(HexMeshing ${VTK_LIBRARIES})

# TODO: Add tests and install targets if needed.
