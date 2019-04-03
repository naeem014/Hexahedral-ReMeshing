#include<vtkGenericDataObjectReader.h>
#include<vtkUnstructuredGridReader.h>
#include<vtkUnstructuredGrid.h>
#include<vtkPolyData.h>
#include<vtkSmartPointer.h>
#include<vtkPoints.h>
#include<vtkCellArray.h>


#include<string>
#include<iostream>
#include<fstream>
#include<limits>
#include "HexMeshing.h"
#include "Mesh.h"

using namespace std;

int main (int argc, char *argv[]) {
	string input_filename  = argv[1];
	vtkSmartPointer<vtkGenericDataObjectReader> reader = vtkSmartPointer<vtkGenericDataObjectReader>::New();
	reader->SetFileName (input_filename.c_str());
    reader->Update();
	if (reader->IsFilePolyData()) {
		vtkPolyData* output = reader->GetPolyDataOutput();
	} else if (reader->IsFileUnstructuredGrid()) {
		vtkUnstructuredGrid* output = reader->GetUnstructuredGridOutput();
		vtkPoints* meshPoints = output->GetPoints();
		Mesh mesh;
		int n_vertices = meshPoints->GetNumberOfPoints();
		mesh.vertices.resize(n_vertices);
		double p[3];
		for (vtkIdType i = 0; i < n_vertices; i++) {
			meshPoints->GetPoint(i, p);
			Vertex v(p[0], p[1], p[2], i);
			mesh.vertices.at(i) = v;
		}
		
		int n_cells = output->GetNumberOfCells();
		mesh.cells.resize(n_cells);
		int n = 0;
		for (vtkIdType i = 0; i < n_cells; i++) {
			vtkSmartPointer<vtkIdList> v_ids = vtkSmartPointer<vtkIdList>::New();
			output->GetCellPoints(i, v_ids);
			vector<int> temp_v_ids;	
			for (vtkIdType j = 0; j < v_ids->GetNumberOfIds(); j++) {
				vtkIdType v_id = v_ids->GetId(j);
				mesh.vertices.at(v_id).cells_ids.push_back(i);
				temp_v_ids.push_back(v_id);		
			}
			mesh.addCell(temp_v_ids, i);
			for (int j = 0; j < mesh.cells.at(i).faces.size(); j++) {
				Face f = mesh.cells.at(i).faces.at(j);
				vtkSmartPointer<vtkIdList> face_v_ids = vtkSmartPointer<vtkIdList>::New();
				for (int k = 0; k < f.v_ids.size(); k++) {
					face_v_ids->InsertNextId(f.v_ids.at(k));			
				}
				vtkSmartPointer<vtkIdList> c_ids = vtkSmartPointer<vtkIdList>::New();
				output->GetCellNeighbors(i, face_v_ids, c_ids);
				if (c_ids->GetNumberOfIds() == 0) {
					mesh.cells.at(i).faces.at(j).isSurfaceFace = true;
                    mesh.surface_faces.push_back(mesh.cells.at(i).faces.at(j));
					mesh.cells.at(i).isSurfaceCell = true;
					mesh.surface_cells.push_back(mesh.cells.at(i));
					mesh.setStepSize(f);
                    mesh.computeNormal(f.v_ids[0], f.v_ids[1], f.v_ids[2]);
					mesh.computeNormal(f.v_ids[1], f.v_ids[0], f.v_ids[2]);
					mesh.computeNormal(f.v_ids[2], f.v_ids[0], f.v_ids[1]);
				} else {
					for (vtkIdType c_id = 0; c_id < c_ids->GetNumberOfIds(); c_id++) {
						mesh.cells.at(i).neighbors.push_back(c_ids->GetId(c_id));
					}
				}
			}			
		}
		mesh.setTypes();
		// vector<Edge> c = mesh.setNeighboringCorners();
		vector<vector<Face>> patches = mesh.extractPatches();
		vector<Edge> c = mesh.extractFacesFromPatch(patches[13]);
		// cout << cvs.size() << endl;
		// vector<Edge> c;
		// for (int i = 0; i < cvs.size(); i+=2) {
		// 	c.emplace_back(cvs.at(i), cvs.at(i+1), c.size());
		// }
		// vector<vector<Face>> patches = mesh.extractPatches();
		// vector<Face> c = patches[0];
		// vector<Edge> c;
		// for (int i = 0; i < patches.size(); i++) {
			// vector<Edge> edges = mesh.extractFacesFromPatch(patches[0]);
			// cout << i << endl;
			// c.insert(c.begin(), edges.begin(), edges.end());
		// }
		// for (int j = 0; j < mesh.corner_cells.size(); j++) {
			// vector<Edge> c = mesh.getBoundaryEdges(0);
			cout << c.size() << endl;

			ofstream outputFile;
			// outputFile.open("output" + to_string(j) + ".vtk");
			outputFile.open("output.vtk");
			outputFile << "# vtk DataFile Version 1.0" << endl;
			outputFile << "Unstructured Grid Example" << endl;
			outputFile << "ASCII" << endl;
			outputFile << "\nDATASET UNSTRUCTURED_GRID" << endl;
			outputFile << "POINTS " << mesh.vertices.size() << " float" << endl;
			for (int i = 0; i < mesh.vertices.size(); i++) {
				outputFile << mesh.vertices.at(i).x << " " << mesh.vertices.at(i).y << " " << mesh.vertices.at(i).z << endl;
			}
			
			int ncells = c.size();
			// int num_points = 0;
			// for (int i = 0; i < ncells; i++) {
			// 	num_points += (1 + c.at(i).v_ids.size());
			// }
			outputFile << "CELLS " << ncells << " " << ncells * 3 << endl;
			for (int i = 0; i < ncells; i++) {
				// outputFile << c.at(i).v_ids.size() << " ";
				// for (int j = 0; j < c.at(i).v_ids.size(); j++) {
				// 	outputFile << c.at(i).v_ids.at(j) << " ";
				// }
				// outputFile << endl;
				// // outputFile << "1 " << c.at(0).v_ids.at(i) << endl;
				outputFile << "2 " << c.at(i).v1 << " " << c.at(i).v2 << endl;
			}


			outputFile << "CELL_TYPES " << ncells << endl;
			for (int i = 0; i < ncells; i++) {
				outputFile << "3" << endl;
			}
		// }
	}
  	return EXIT_SUCCESS;
}
