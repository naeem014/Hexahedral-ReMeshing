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
#include "vtkImageData.h"
#include "vtkCoordinate.h"
#include "vtkSphereSource.h"
#include "vtkButtonWidget.h"
#include "vtkTexturedButtonRepresentation2D.h"

#include<string>
#include<iostream>
#include<fstream>
#include<limits>
#include "HexMeshing.h"
#include "Mesh.h"

using namespace std;

vector<unsigned char> readPNG(const char* filename) {
	vector<unsigned char> image;
	unsigned width, height;
	return image;
}

void CreateImage(vtkSmartPointer<vtkImageData> image,
                 unsigned char* color1,
                 unsigned char* color2)
{
  // Specify the size of the image data
  image->SetDimensions(10, 10, 1);
  image->AllocateScalars(VTK_UNSIGNED_CHAR, 3);

  int* dims = image->GetDimensions();

  // Fill the image with
  for (int y = 0; y < dims[1]; y++)
  {
    for (int x = 0; x < dims[0]; x++)
    {
      unsigned char* pixel =
        static_cast<unsigned char*>(image->GetScalarPointer(x, y, 0));
      if (x < 5)
      {
        pixel[0] = color1[0];
        pixel[1] = color1[1];
        pixel[2] = color1[2];
      }
      else
      {
        pixel[0] = color2[0];
        pixel[1] = color2[1];
        pixel[2] = color2[2];
      }
    }
  }
}

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
		renWin->FullScreenOn();
		// renWin->SetCurrentCursor(VTK_CURSOR_CROSSHAIR);

		vtkSmartPointer<vtkImageData> image1 = vtkSmartPointer<vtkImageData>::New();
		vtkSmartPointer<vtkImageData> image2 =	vtkSmartPointer<vtkImageData>::New();
		unsigned char banana[3] = { 227, 207, 87 };
		unsigned char tomato[3] = { 255, 99, 71 };
		CreateImage(image1, banana, tomato);
		CreateImage(image2, tomato, banana);

		vtkRenderWindowInteractor *iren = vtkRenderWindowInteractor::New();

		vtkSmartPointer<vtkTexturedButtonRepresentation2D> buttonRepresentation =
			vtkSmartPointer<vtkTexturedButtonRepresentation2D>::New();
		buttonRepresentation->SetNumberOfStates(2);
		buttonRepresentation->SetButtonTexture(0, image1);
		buttonRepresentation->SetButtonTexture(1, image2);

		vtkSmartPointer<vtkButtonWidget> buttonWidget =
			vtkSmartPointer<vtkButtonWidget>::New();
		buttonWidget->SetInteractor(iren);
		buttonWidget->SetRepresentation(buttonRepresentation);

    	iren->SetRenderWindow(renWin);

		ren1->SetBackground(0.8,0.3,0.3);

		renWin->Render();

		vtkSmartPointer<vtkCoordinate> upperRight =
			vtkSmartPointer<vtkCoordinate>::New();
		upperRight->SetCoordinateSystemToNormalizedDisplay();
		upperRight->SetValue(1.0, 1.0);

		double bds[6];
		double sz = 50.0;
		bds[0] = upperRight->GetComputedDisplayValue(ren1)[0] - sz;
		bds[1] = bds[0] + sz;
		bds[2] = upperRight->GetComputedDisplayValue(ren1)[1] - sz;
		bds[3] = bds[2] + sz;
		bds[4] = bds[5] = 0.0;

		// Scale to 1, default is .5
		buttonRepresentation->SetPlaceFactor(1);
		buttonRepresentation->PlaceWidget(bds);

		buttonWidget->On();


		iren->Initialize();
    	iren->Start();
		iren->StartPinchEvent();*/
		vtkPoints* meshPoints = output->GetPoints();
		Mesh mesh;
		int n_vertices = meshPoints->GetNumberOfPoints();
		mesh.vertices.resize(n_vertices);
		double p[3];
		double scaleUpFactor = 1;
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
		// mesh.setStepSizes();
		// vector<vector<Face>> patches = mesh.extractPatches();
		// vector<Edge> c = mesh.extractFacesFromPatch(patches[8]);
		// cout << cvs.size() << endl;
		// vector<Edge> c = mesh.corner_edges;
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
			// vector<int> c;
			vector<Edge> c = mesh.getBoundaryEdges();
			
			for (int i = 0; i < mesh.corner_edges.size(); i++) {
				Vertex v1 = mesh.vertices.at(mesh.corner_edges.at(i).v1);
				int id1 = mesh.vertices.size();
				mesh.vertices.emplace_back(v1.u, v1.v, v1.w, id1);
				Vertex v2 = mesh.vertices.at(mesh.corner_edges.at(i).v2);
				int id2 = mesh.vertices.size();
				mesh.vertices.emplace_back(v2.u, v2.v, v2.w, id2);
				c.emplace_back(id1, id2, c.size());
				
			// 	for (int j = 0; j < v1.neighboring_corners.size(); j++) {
			// 		Vertex& v2 = mesh.vertices.at(v1.neighboring_corners.at(j));
			// 		int id2 = mesh.vertices.size();
			// 		mesh.vertices.emplace_back(v2.u, v2.v, v2.w, id2);
			// 		c.emplace_back(id1, id2, c.size());
			// 	}
			// 	vector<vector<double>> mask_points = mesh.getMaskCoords(v);
			// 	for (int j = 0; j < mask_points.size(); j++) {
			// 		vector<double> p = mask_points.at(j);
			// cout << "x: " << v.x << " y: " << v.y << " z: " << v.z << endl;
			// cout << "u: " << v.u << " v: " << v.v << " w: " << v.w << endl;
			// cout << "=======================================================" << endl;
 				
			}
			// 	break;
			// }
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
			outputFile << "CELLS " << ncells << " " << ncells * 3 << endl;
			for (int i = 0; i < ncells; i++) {
				// outputFile << c.at(i).v_ids.size() << " ";
				// for (int j = 0; j < c.at(i).v_ids.size(); j++) {
				// 	outputFile << c.at(i).v_ids.at(j) << " ";
				// }
				// outputFile << endl;
				// outputFile << "1 " << c.at(i) << endl;
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
