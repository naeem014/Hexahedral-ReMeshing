### Hexahedral-ReMeshing

Remeshing Tetrahedral Polycube models to Hexahedral mesh. This project takes a vtk file of a tetrahedral polycube model and created a mesh data structure. The mesh data structure implements functions to extract boundary surfaces of input mesh, extract boundary edges of input mesh, divides each surface into quads, translates corner points to regular integer grid and volume division.

## Project setup
- Clone the project to local machine.
- Enter the project directory and make new directory e.g mkdir build
- Open CMake GUI which can be downloaded from here: https://cmake.org/download/
- Click Browse source and select the project directory
- Click Browse build and select the build folder inside project directory
- Click Configure
- Click Generate
- Click Open Project (The project will be opened in Visual Studio. The latest version of Visual Studio is recommended.)
- Build the project in Visual Studio
- Open terminal and navigate to build/Debug inside the project directory
- Run ./HexMeshing ../../Examples/kitty/polycube.tet.vtk
- The above command outputs two files: surface_division.vtk and volume_division.vtk
- Paraview can be used to open the vtk files to visualize the results.
