#include<vtkGenericDataObjectReader.h>
#include<vtkUnstructuredGridReader.h>
#include<vtkUnstructuredGrid.h>
#include<vtkPolyData.h>
#include<vtkSmartPointer.h>
#include<vtkPoints.h>
#include<vtkCellArray.h>

#include "vtkDataSetMapper.h"
#include "vtkProperty.h"
#include "vtkVolumeProperty.h"
#include "vtkColorTransferFunction.h"
#include "vtkPiecewiseFunction.h"
#include "vtkSphereSource.h"
#include "vtkPolyDataMapper.h"
#include "vtkActor.h"
#include "vtkRenderWindow.h"
#include "vtkRenderer.h"
#include "vtkRenderWindowInteractor.h"
#include "vtkUnstructuredGridVolumeRayCastMapper.h"
#include "vtkDataSetTriangleFilter.h"
#include "vtkUnstructuredGridVolumeMapper.h"
#include "vtkProjectedTetrahedraMapper.h"
#include "vtkFloatArray.h"
#include "vtkPointData.h"

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
		

		/*vtkDataSetMapper *meshMapper = vtkDataSetMapper::New();
		meshMapper->SetInputData(output);

		vtkSmartPointer<vtkActor> meshActor = vtkSmartPointer<vtkActor>::New();
		meshActor->SetMapper(meshMapper);

		vtkRenderer *ren1 = vtkRenderer::New();
		vtkRenderWindow *renWin = vtkRenderWindow::New();
		renWin->AddRenderer(ren1);
		ren1->AddActor(meshActor);

		vtkRenderWindowInteractor *iren = vtkRenderWindowInteractor::New();
    	iren->SetRenderWindow(renWin);

		ren1->SetBackground(1,1,1);

		renWin->Render();

		iren->Initialize();
    	iren->Start();*/

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
		mesh.setNeighboringCorners();
		// vector<vector<Face>> patches = mesh.extractPatches();
		// vector<Edge> c = mesh.extractFacesFromPatch(patches[8]);
		// cout << cvs.size() << endl;
		// vector<Edge> c;
		// for (int i = 0; i < mesh.corner_vertices.size(); i++) {
		// 	Vertex& v = mesh.vertices.at(mesh.corner_vertices.at(i));
		// 	cout << v.neighboring_corners.size() << endl;
		// 	for (int j = 0; j < v.neighboring_corners.size(); j++) {
		// 		c.emplace_back(v.id, mesh.vertices.at(v.neighboring_corners.at(j)).id, c.size());
		// 	}
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
			// cout << mesh.stepSize << endl;
			// vector<Edge> c = mesh.getBoundaryEdges();
			// cout << "Did I actually arrive here?" << endl;
			vector<int> c;
			for (int i = 0; i < mesh.corner_vertices.size(); i++) {
				Vertex& v = mesh.vertices.at(mesh.corner_vertices.at(i));
				vector<vector<double>> mask_points = mesh.getMaskCoords(v);
				for (int j = 0; j < mask_points.size(); j++) {
					vector<double> p = mask_points.at(j);
					int id = mesh.vertices.size();
					mesh.vertices.emplace_back(p[0], p[1], p[2], id);
					c.push_back(id);
				}
				break;
			}
			// vector<Face> c = mesh.mesh_faces;
			// int num_points = 0;
			// for (int i = 0; i < c.size(); i++) {
			// 	for (int j = 0; j < c.at(i).v_ids.size(); j++) {
			// 		num_points += 1;
			// 	}
			// 	num_points += 1;
			// }
			// cout << num_points << endl;
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
			outputFile << "CELLS " << ncells << " " << ncells * 2 << endl;
			for (int i = 0; i < ncells; i++) {
				// outputFile << c.at(i).v_ids.size() << " ";
				// for (int j = 0; j < c.at(i).v_ids.size(); j++) {
				// 	outputFile << c.at(i).v_ids.at(j) << " ";
				// }
				// outputFile << endl;
				outputFile << "1 " << c.at(i) << endl;
				// outputFile << "2 " << c.at(i).v1 << " " << c.at(i).v2 << endl;
			}


			outputFile << "CELL_TYPES " << ncells << endl;
			for (int i = 0; i < ncells; i++) {
				outputFile << "1" << endl;
			}
		// }
	// }
  	return EXIT_SUCCESS;
}
