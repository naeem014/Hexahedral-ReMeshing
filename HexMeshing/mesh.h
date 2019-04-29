#ifndef __MESH_H__
#define __MESH_H__

#include<string>
#include<vector>
#include <algorithm>
#include "glm/glm.hpp"

using namespace std;

class Vertex {
    public:
        double x;
        double y;
        double z;
        int id;
        int normal_x = 0;
        int normal_y = 0;
        int normal_z = 0;
        double u, v, w;
        vector<double> mask_coords;
        int normal_sum = 0;
        int visitCount = 0;
        vector<int> neighbors;
        vector<int> neighboring_corners;
        vector<int> cells_ids;
        vector<int> edge_ids;
        vector<int> corner_neighbors;
        bool isSurfaceVertex = false;
        bool isBoundaryVertex = false;
        bool isCornerVertex = false;
        bool isVisited = false;
        Vertex() {}
        Vertex(double x_, double y_, double z_, int id_) {
            x = x_;
            y = y_;
            z = z_;
            id = id_;
        }
};

class Edge {
    public:
        int v1, v2;
        int id;
        bool isBoundaryEdge = false;
        bool isCornerEdge = false;
        bool isVisited = false;
        vector<int> neighbors;
        string type;
        Edge() {}
        Edge(int id1, int id2, int id_) {
            v1 = id1;
            v2 = id2;
            id = id_;
        }
};

class Face {
    public:
        vector<int> v_ids;
        vector<int> neighbors;
        bool isBoundaryFace = false;
        bool isSurfaceFace = false;
        bool isCornerFace = false;
        bool isVisited = false;
        Face() {}
        Face(int v_id0, int v_id1, int v_id2) {
            v_ids.push_back(v_id0);
            v_ids.push_back(v_id1);
            v_ids.push_back(v_id2);
        }

        Face(vector<int> ids) {
            for (int i = 0; i < ids.size(); i++) {
                v_ids.push_back(ids[i]);
            }
        }
};

class Cell {
    public:
        vector<Face> faces;
        vector<int> v_ids;
        vector<int> neighbors;
        int id;
        bool isVisited = false;
        bool isSurfaceCell = false;
        bool isBoundaryCell = false;
        bool isCornerCell = false;
};

class Mesh {
    public:
        double PI = 3.14159265;
        vector<Vertex> vertices;
        vector<Cell> cells, surface_cells, boundary_cells, corner_cells;
        vector<Face> faces, surface_faces, boundary_faces, corner_faces, mesh_faces;
        vector<int> surface_vertices, boundary_vertices, corner_vertices;
        vector<Edge> corner_edges;
        double stepSize = numeric_limits<double>::infinity();
        double stepSizeX = numeric_limits<double>::infinity();
        double stepSizeY = numeric_limits<double>::infinity();
        double stepSizeZ = numeric_limits<double>::infinity();
        void addCell(vector<int> v_ids, int loc);
        void setTypes();
        void computeNormal(int v_id0, int v_id1, int v_id2);
        vector<Face> extractFaces();
        vector<vector<Face>> extractPatches();
        void setNeighboringCorners();
        void setCornerNeighbors(vector<Face> patch);
        vector<Edge> extractFacesFromPatch(vector<Face> patch);
        vector<double> getPlaneNormal(Face& f);
        bool isPlaneSame(Face& f, vector<double> plane_normal);
        bool edgeInNeighbors(Edge e, vector<Face> patch, int id);
        bool containsEdges(vector<int> v_ids1, vector<int> v_ids2);
        bool containsVertices(vector<int> v_ids1, vector<int> v_ids2);
        void setStepSize(Face& f);
        void setStepSizes();
        void setCornersUVWcoords();
        vector<double> getCrossProduct(vector<double> a, vector<double> b);
        double getDotProduct(vector<double> a, vector<double> b);
        double getVectorLength(vector<double> v);
        vector<double> getDirectionVector(Vertex v1, Vertex v2);
        double doubleDifference(double a, double b);
        vector<double> normalizeVector(vector<double> v);
        vector<double> getBarycentricCoordinates(Vertex p, Face f);
        vector<double> getBarycentricCoordinates(Vertex p, int cell_id);
        vector<Edge> getBoundaryEdges();
        int getCornerIndex(Face& f);
        vector<Edge> traceLineFromCorner(int corner_id, int c_id, vector<double> direction);
        vector<double> rotateVector(vector<double> v, double angle, int axis);
        vector<vector<double>> getMaskCoords(Vertex& v);
        bool isPointInsideVolume(vector<double> p);
        bool isPointOnPlane(vector<double> p, Face& plane, glm::vec3 normal, vector<glm::vec3> directions);
        bool isPointOutside(vector<double> p, Vertex& v);
        void alignNeighbors(Vertex& v, int direction);
};

#endif