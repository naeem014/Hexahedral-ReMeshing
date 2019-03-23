# CMakeList.txt : CMake project for HexMeshing, include source and define
# project specific logic here.
#
cmake_minimum_required (VERSION 2.8)

find_package(VTK REQUIRED)
include(${VTK_USE_FILE})
include_directories(
	${VTK_INCLUDE_DIRS},
	${PROJECT_SOURCE_DIR}/glm/glm
)
# Add source to this project's executable.
set(EXT_LIBS ${VTK_LIBRARIES})
set(SRC_HEX_MESHING
	mesh.cpp
	mesh.h
)
add_executable (HexMeshing ${SRC_HEX_MESHING} "HexMeshing.cpp" "HexMeshing.h")

target_link_libraries(HexMeshing ${EXT_LIBS})

# TODO: Add tests and install targets if needed.
