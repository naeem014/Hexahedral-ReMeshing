#include<iostream>
#include<math.h>
#include<algorithm>
#include<limits>
#include "glm/glm.hpp"

#include "Mesh.h"

#define PI 3.141459265

using namespace std;

void print_vector(string message, vector<double> v);

void Mesh::addCell(vector<int> v_ids, int loc) {
    Cell c;
    for (int i = 0; i < v_ids.size(); i++) {
        c.v_ids.push_back(v_ids.at(i));
        for (int j = 0; j < v_ids.size(); j++) {
            if (i != j && find(vertices.at(v_ids[i]).neighbors.begin(), vertices.at(v_ids[i]).neighbors.end(), vertices.at(v_ids[j]).id) == vertices.at(v_ids[i]).neighbors.end()) {
                vertices.at(v_ids[i]).neighbors.push_back(vertices.at(v_ids[j]).id);
            }
        }
    }
    Face f1(v_ids[0], v_ids[1], v_ids[2]);
    c.faces.push_back(f1);
    Face f2(v_ids[0], v_ids[1], v_ids[3]);
    c.faces.push_back(f2);
    Face f3(v_ids[0], v_ids[2], v_ids[3]);
    c.faces.push_back(f3);
    Face f4(v_ids[1], v_ids[2], v_ids[3]);
    c.faces.push_back(f4);
    c.id = loc;
    cells.at(loc) = c;
}

void Mesh::setTypes() {
    for (int i = 0; i < surface_cells.size(); i++) {
        Cell& c = cells.at(surface_cells.at(i).id);
        for (int j = 0; j < c.faces.size(); j++) {
            if (c.faces.at(j).isSurfaceFace) {
                Face& f = c.faces.at(j);
                for (int l = 0; l < f.v_ids.size(); l++) {
                    vertices.at(f.v_ids[l]).isSurfaceVertex = true;
                    surface_vertices.push_back(f.v_ids[l]);
                    if (vertices.at(f.v_ids[l]).normal_sum == 2) {
                        vertices.at(f.v_ids[l]).isBoundaryVertex = true;
                        boundary_vertices.push_back(f.v_ids[l]);
                        cells.at(surface_cells.at(i).id).faces.at(j).isBoundaryFace = true;
                        boundary_faces.push_back(cells.at(surface_cells.at(i).id).faces.at(j));
                        cells.at(surface_cells.at(i).id).isBoundaryCell = true;
                        boundary_cells.push_back(cells.at(surface_cells.at(i).id));
                    } else if (vertices.at(f.v_ids[l]).normal_sum == 3) {
                        vertices.at(f.v_ids[l]).isBoundaryVertex = true;
                        boundary_vertices.push_back(f.v_ids[l]);
                        cells.at(surface_cells.at(i).id).faces.at(j).isBoundaryFace = true;
                        boundary_faces.push_back(cells.at(surface_cells.at(i).id).faces.at(j));
                        cells.at(surface_cells.at(i).id).isBoundaryCell = true;
                        boundary_cells.push_back(cells.at(surface_cells.at(i).id));
                        vertices.at(f.v_ids[l]).isCornerVertex = true;
                        corner_vertices.push_back(f.v_ids[l]);
                        cells.at(surface_cells.at(i).id).faces.at(j).isCornerFace = true;
                        corner_faces.push_back(cells.at(surface_cells.at(i).id).faces.at(j));
                        cells.at(surface_cells.at(i).id).isCornerCell = true;
                        corner_cells.push_back(cells.at(surface_cells.at(i).id));
                    }
                }
            }
        }
    }
    for (int i = 0; i < corner_vertices.size(); i++) {
        if (count(corner_vertices.begin(), corner_vertices.end(), corner_vertices.at(i)) > 1) {
            corner_vertices.erase(corner_vertices.begin() + i);
            i = 0;
        }
    }
}

vector<vector<Face>> Mesh::extractPatches() {
    vector<vector<Face>> patches;
    vector<Cell> queue = boundary_cells;
    vector<double> plane_normal;
    int src_id = -1;
    while (!queue.empty()) {
        Cell& c = cells.at(queue.back().id);
        bool faceFound = false;
        int f_id = -1;
        for (int i = 0; i < c.faces.size(); i++) {
            if (c.faces.at(i).isSurfaceFace && !c.faces.at(i).isVisited && isPlaneSame(c.faces.at(i), plane_normal)) {
                faceFound = true;
                f_id = i;
                c.faces.at(i).isVisited = true;
                patches.back().push_back(c.faces.at(i));
                plane_normal = getPlaneNormal(c.faces.at(i));
                break;
            }
        }
        if (!faceFound) {
            for (int i = 0; i < c.faces.size(); i++) {
                if (c.faces.at(i).isSurfaceFace && !c.faces.at(i).isVisited) {
                    plane_normal = getPlaneNormal(c.faces.at(i));
                    faceFound = true;
                    f_id = i;
                    c.faces.at(i).isVisited = true;
                    patches.emplace_back();
                    patches.back().push_back(c.faces.at(i));
                    src_id = c.id;
                    break;
                }
            }
            if (!faceFound) {
                queue.pop_back();
                if (src_id == c.id) {
                    plane_normal.clear();
                }
                continue;
            }
        }
        for (int i = 0; i < c.faces.at(f_id).v_ids.size(); i++) {
            Vertex& v = vertices.at(c.faces.at(f_id).v_ids.at(i));
            for (int j = 0; j < v.cells_ids.size(); j++) {
                Cell& c_neighbor = cells.at(v.cells_ids.at(j));
                if (c_neighbor.isSurfaceCell && c_neighbor.id != c.id) {
                    for (int k = 0; k < c.faces.size(); k++) {
                        Face& f_neighbor = c_neighbor.faces.at(k);
                        if (f_neighbor.isSurfaceFace && !f_neighbor.isVisited && isPlaneSame(f_neighbor, plane_normal)) {
                            queue.push_back(c_neighbor);
                            break;
                        }
                    }
                }
            }
        }
    }
    return patches;
}

bool Mesh::isPlaneSame(Face& f, vector<double> plane_normal) {
    if (plane_normal.empty()) {
        return false;
    } else {
        vector<double> normal = getPlaneNormal(f);
        double dotProduct = getDotProduct(normal, plane_normal);
        double angle = getDotProduct(normal, plane_normal) / (getVectorLength(normal) * getVectorLength(plane_normal));
        if (angle > 0.9) {
            return true;
        }
    }
    return false;
}

vector<double> Mesh::getPlaneNormal(Face& f) {
    
    Vertex v0 = vertices.at(f.v_ids[0]);
    Vertex v1 = vertices.at(f.v_ids[1]);
    Vertex v2 = vertices.at(f.v_ids[2]);

    double a[] = {v1.x - v0.x, v1.y - v0.y, v1.z - v0.z};
    double b[] = {v2.x - v0.x, v2.y - v0.y, v2.z - v0.z};

    double normal_x = (a[1] * b[2]) - (a[2] * b[1]);
    double normal_y = (a[2] * b[0]) - (a[0] * b[2]);
    double normal_z = (a[0] * b[1]) - (a[1] * b[0]);

    double length = sqrt((normal_x * normal_x) + (normal_y * normal_y) + (normal_z * normal_z));

    vector<double> normal = {abs(normal_x / length), abs(normal_y / length), abs(normal_z / length)};
    return normal;
}

vector<Edge> Mesh::getBoundaryEdges() {
    vector<Edge> boundary_edges;
    vector<glm::vec3> directions = {
        glm::vec3(1, 0, 0),
        glm::vec3(0, 1, 0),
        glm::vec3(0, 0, 1),
    };
    vector<Edge> edges;
    double step;
    for (int i = 0; i < corner_edges.size(); i++) {
        // boundary_edges.push_back(corner_edges.at(i));
        Vertex v1 = vertices.at(corner_edges.at(i).v1);
        Vertex v2 = vertices.at(corner_edges.at(i).v2);
        glm::vec3 d1 = glm::normalize(glm::vec3(v2.u - v1.u, v2.v - v1.v, v2.w - v1.w));
        glm::vec3 d2 = glm::normalize(glm::vec3(v1.u - v2.u, v1.v - v2.v, v1.w - v2.w));
        
        if (round(abs(glm::dot(d1, directions[0]))) == 1) {
            step = stepSizeX;
        } else if (round(abs(glm::dot(glm::normalize(d1), directions[1]))) == 1) {
            step = stepSizeY;
        } else if (round(abs(glm::dot(glm::normalize(d1), directions[2]))) == 1) {
            step = stepSizeZ;
        }
        vector<double> direction1 = {d1[0], d1[1], d1[2]};
        vector<double> direction2 = {d2[0], d2[1], d2[2]};
        
        edges = traceLineFromCorner(v1.id, direction2, step);
        boundary_edges.insert(boundary_edges.begin(), edges.begin(), edges.end());
        edges.clear();
        edges = traceLineFromCorner(v2.id, direction1, step);
        boundary_edges.insert(boundary_edges.begin(), edges.begin(), edges.end());
        edges.clear();
    }
    
    return boundary_edges;
}

vector<Edge> Mesh::traceLineFromCorner(int v_id, vector<double> direction, double step) {
    vector<Edge> boundary_edges;
    Vertex& corner = vertices.at(v_id);
    Cell& c = cells.at(corner.cells_ids.at(0));

    glm::vec3 initial(corner.u, corner.v, corner.w);
    glm::vec3 d(step * direction[0], step * direction[1], step * direction[2]);

    int c_id = -1;
    while (true) {
        vertices.emplace_back(initial[0], initial[1], initial[2], vertices.size());
        Vertex start = vertices.at(vertices.size() - 1);
        glm::vec3 new_end = initial + d;
        vertices.emplace_back(new_end[0], new_end[1], new_end[2], vertices.size());
        Vertex end = vertices.at(vertices.size() - 1);
        boundary_edges.emplace_back(start.id, end.id, boundary_edges.size());

        if (c_id == -1) {
            bool foundCell = false;
            for (int i = 0; i < c.v_ids.size(); i++) {
                Vertex& c_v = vertices.at(c.v_ids.at(i));
                for (int j = 0; j < c_v.cells_ids.size(); j++) {
                    if (c_v.cells_ids.at(j) == c_id) {
                        continue;
                    }
                    vector<double> b_coords = getBarycentricCoordinates(end, c_v.cells_ids.at(j));
                    if (b_coords[0] >= 0 && b_coords[0] <= 1 && b_coords[1] >= 0 && b_coords[1] <= 1
                        && b_coords[2] >= 0 && b_coords[2] <= 1 && b_coords[3] >= 0 && b_coords[3] <= 1) {
                        initial = new_end;
                        c_id = c_v.cells_ids.at(j);
                        c = cells.at(c_id);
                        foundCell = true;
                        break;
                    }
                }
                if (foundCell) {
                    break;
                }
            }
            if (foundCell) {
                continue;
            } else {
                break;
            }
        }
        vector<double> b_coords = getBarycentricCoordinates(end, c_id);
        if (b_coords[0] >= 0 && b_coords[0] <= 1 && b_coords[1] >= 0 && b_coords[1] <= 1
            && b_coords[2] >= 0 && b_coords[2] <= 1 && b_coords[3] >= 0 && b_coords[3] <= 1) {
            // cout << b_coords[0] << " " << b_coords[1] << " " << b_coords[2] << " " << b_coords[3] << endl;
            initial = new_end;
            continue;
        } else {
            // cout << "==========================================" << endl;
            c_id = -1;
            continue;
        }
        break;
    }
    return boundary_edges;
}

vector<vector<double>> Mesh::getMaskCoords(Vertex& v) {
    vector<vector<double>> points;
    points.push_back({v.x + (stepSize / 2), v.y + (stepSize / 2), v.z + (stepSize / 2)}); 
    points.push_back({v.x + (stepSize / 2), v.y + (stepSize / 2), v.z - (stepSize / 2)}); 
    points.push_back({v.x + (stepSize / 2), v.y - (stepSize / 2), v.z + (stepSize / 2)}); 
    points.push_back({v.x + (stepSize / 2), v.y - (stepSize / 2), v.z - (stepSize / 2)}); 
    points.push_back({v.x - (stepSize / 2), v.y + (stepSize / 2), v.z + (stepSize / 2)}); 
    points.push_back({v.x - (stepSize / 2), v.y + (stepSize / 2), v.z - (stepSize / 2)}); 
    points.push_back({v.x - (stepSize / 2), v.y - (stepSize / 2), v.z + (stepSize / 2)});
    points.push_back({v.x - (stepSize / 2), v.y - (stepSize / 2), v.z - (stepSize / 2)});
    // points.push_back({v.x + (stepSize / 2), v.y, v.z}); 
    // points.push_back({v.x - (stepSize / 2), v.y, v.z}); 
    // points.push_back({v.x, v.y + (stepSize / 2), v.z}); 
    // points.push_back({v.x, v.y - (stepSize / 2), v.z}); 
    // points.push_back({v.x, v.y, v.z + (stepSize / 2)}); 
    // points.push_back({v.x, v.y, v.z - (stepSize / 2)}); 
    vector<vector<double>> new_points;
    int index = 0;
    for (int i = 0; i < points.size(); i++) {
        if (isPointInsideVolume(points.at(i))) {
            // index = i;
            // break;
            new_points.push_back(points.at(i));
        }
    }
    // return points.at(index);
    return new_points;
}

bool Mesh::isPointInsideVolume(vector<double> p) {
    vector<glm::vec3> directions = {
        glm::vec3(1, 0, 0),
        glm::vec3(-1, 0, 0),
        glm::vec3(0, 1, 0),
        glm::vec3(0, -1, 0),
        glm::vec3(0, 0, 1),
        glm::vec3(0, 0, -1)
    };
    for (int i = 0; i < directions.size(); i++) {
        int count = 0;
        glm::vec3 d = glm::vec3(directions.at(i)[0], directions.at(i)[1], directions.at(i)[2]);
        // cout << "Direction at " << i << ": " << d[0] << " " << d[1] << " " << d[2] << endl;
        for (int j = 0; j < mesh_faces.size(); j++) {
            // cout << "Face at " << j << endl; 
            bool doesLineIntersectPlane = false;
            Face& plane = mesh_faces.at(j);
            Vertex& p1 = vertices.at(plane.v_ids.at(1));
            Vertex& p2 = vertices.at(plane.v_ids.at(2));
            Vertex& p3 = vertices.at(plane.v_ids.at(0));
            glm::vec3 normal = glm::cross(glm::vec3(p2.x - p1.x, p2.y - p1.y, p2.z - p1.z), glm::vec3(p3.x - p1.x, p3.y - p1.y, p3.z - p1.z));
            glm::vec3 PN = glm::vec3(p[0] - p1.x, p[1] - p1.y, p[2] - p1.z);
            normal = glm::normalize(normal);
            double dotPN = glm::dot(normal, PN);
            // cout << "PN: " << PN[0] << " " << PN[1] << " " << PN[2] << endl;
            // cout << "normal: " << n[0] << " " << n[1] << " " << n[2] << endl;
            // cout << "Plane P dot: " << dotPN << endl;
            if (dotPN == 0) {
                doesLineIntersectPlane = isPointOnPlane(p, plane, normal, directions);
            } else {
                double A = normal[0];
                double B = normal[1];
                double C = normal[2];
                double D = - (A * p1.x) - (B * p1.y) - (C * p1.z);
                double t = - ((A * p[0]) + (B * p[1]) + (C * p[2]) + D) / ((A * d[0]) + (B * d[1]) + (C * d[2]));
                vector<double> new_p = {p[0] + (t * d[0]), p[1] + (t * d[1]), p[2] + (t * d[2])};
                doesLineIntersectPlane = isPointOnPlane(new_p, plane, normal, directions);
            }
            if (doesLineIntersectPlane) {
                count += 1;
            }
        }
        if (count % 2 == 0) {
            return false;
        }
        // cout << "==============================================================================" << endl;
    }
    return true;
}

bool Mesh::isPointOnPlane(vector<double> p, Face& plane, glm::vec3 normal, vector<glm::vec3> directions) {
    for (int i = 0; i < directions.size(); i++) {
        int count = 0;
        if (abs(glm::dot(directions.at(i), glm::normalize(normal))) == 1) {
            continue;
        }
        for (int j = 0; j < plane.v_ids.size(); j++) {
            Vertex& a = vertices.at(plane.v_ids.at(j));
            Vertex& b = vertices.at(plane.v_ids.at((j + 1) % plane.v_ids.size()));
            glm::vec3 p1 = glm::vec3(a.x, a.y, a.z);
            glm::vec3 v1 = glm::vec3(b.x - a.x, b.y - a.y, b.z - a.z);
            glm::vec3 p2 = glm::vec3(p[0], p[1], p[2]);
            glm::vec3 v2 = directions.at(i);
            glm::vec3 p3 = glm::vec3(b.x, b.y, b.z);

            glm::vec3 p2_p1 = p2 - p1;
            glm::vec3 v1Xv2 = glm::cross(v1, v2);
            glm::vec3 p2_p1Xv2 = glm::cross(p2_p1, v2);
            
            double A = glm::length(p2_p1Xv2) / glm::length(v1Xv2);
            // cout << glm::dot(p2_p1Xv2, v1Xv2) << endl;
            if (glm::dot(p2_p1Xv2, v1Xv2) < 0) {
                A = -A;
                continue;
            }
            glm::vec3 intersection = p1 + glm::vec3(A * v1[0], A * v1[1], A * v1[2]);
            if ((glm::length(p1 - intersection) / glm::length(v1)) <= 1 && (glm::length(p3 - intersection) / glm::length(v1)) <= 1) {
                count += 1;
            }
        }
        // cout << count << endl;
        if (count % 2 == 0) {
            return false;
        }
    }
    return true;
}

int Mesh::getCornerIndex(Face& f) {
    int index = 0;
    for (int i = 0; i < f.v_ids.size(); i++) {
        if (vertices.at(f.v_ids[i]).isCornerVertex) {
            index = i;
        }
    }
    return index;
}

vector<Edge> Mesh::traceLineFromCorner(int corner_id, int c_id, vector<double> direction) {
    vector<Edge> boundary_edges;
    Cell& c = cells.at(c_id);
    Vertex& corner = vertices.at(corner_id);
    double step = 1;
    // step = stepSize;
    glm::vec3 initial(corner.mask_coords[0], corner.mask_coords[1], corner.mask_coords[2]);
    glm::vec3 d(step * direction[0], step * direction[1], step * direction[2]);
    // cout << d[0] << " " << d[1] << " " << d[2] << endl;
    while (true) {
        int id = vertices.size();
        Vertex start(initial[0], initial[1], initial[2], id);
        vertices.push_back(start);
        id = vertices.size();
        glm::vec3 new_end = initial + d;
        Vertex end(new_end[0], new_end[1], new_end[2], id);
        vertices.push_back(end);
        boundary_edges.emplace_back(start.id, end.id, boundary_edges.size());
        vector<double> b_coords = getBarycentricCoordinates(end, c_id);
        // cout << "C_ID: " << c_id << endl;
        // cout << "a: " << b_coords[0] << " b: " << b_coords[1] << " c: " << b_coords[2] << endl;// << " d: " << b_coords[3] << endl;
        // cout << b_coords[0] + b_coords[1] + b_coords[2] << endl;
        // cout << c_id << endl;
        // cout << "======================================================================================" << endl;
        if (b_coords[0] >= 0 && b_coords[0] <= 1 && b_coords[1] >= 0 && b_coords[1] <= 1 && b_coords[2] >= 0 && b_coords[2] <= 1) {// && b_coords[3] >= 0) {
            initial = new_end;
            continue;
        } else {
            bool foundCell = false;
            for (int i = 0; i < c.v_ids.size(); i++) {
                Vertex& c_v = vertices.at(c.v_ids.at(i));
                for (int j = 0; j < c_v.cells_ids.size(); j++) {
                    if (c_v.cells_ids.at(j) == c_id) {
                        continue;
                    }
                    b_coords = getBarycentricCoordinates(end, c_v.cells_ids.at(j));
                    if (b_coords[0] >= 0 && b_coords[0] <= 1 && b_coords[1] >= 0 && b_coords[1] <= 1 && b_coords[2] >= 0 && b_coords[2] <= 1) {// && b_coords[3] >= 0) {
                        // cout << "a: " << b_coords[0] << " b: " << b_coords[1] << " c: " << b_coords[2] << endl;// << " d: " << b_coords[3] << endl;
                        // cout << "||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||" << endl;
                        initial = new_end;
                        c_id = c_v.cells_ids.at(j);
                        c = cells.at(c_id);
                        foundCell = true;
                        break;
                    }
                }
                if (foundCell) {
                    break;
                }
            }
            if (foundCell) {
                continue;
            }
        }
        break;
    }
    return boundary_edges;
}

vector<double> Mesh::getBarycentricCoordinates(Vertex p, int cell_id) {
    Cell& cell = cells.at(cell_id);
    
    Vertex a = vertices.at(cell.v_ids[0]);
    Vertex b = vertices.at(cell.v_ids[1]);
    Vertex c = vertices.at(cell.v_ids[2]);
    Vertex d = vertices.at(cell.v_ids[3]);

    glm::vec3 vap = glm::vec3(p.x - a.x, p.y - a.y, p.z - a.z);
    glm::vec3 vbp = glm::vec3(p.x - b.x, p.y - b.y, p.z - b.z);
    glm::vec3 vab = glm::vec3(b.x - a.x, b.y - a.y, b.z - a.z);
    glm::vec3 vac = glm::vec3(c.x - a.x, c.y - a.y, c.z - a.z);
    glm::vec3 vad = glm::vec3(d.x - a.x, d.y - a.y, d.z - a.z);
    glm::vec3 vbc = glm::vec3(c.x - b.x, c.y - b.y, c.z - b.z);
    glm::vec3 vbd = glm::vec3(d.x - b.x, d.y - b.y, d.z - b.z);

    double va6 = glm::dot(vbp, glm::cross(vbd, vbc));
    double vb6 = glm::dot(vap, glm::cross(vac, vad));
    double vc6 = glm::dot(vap, glm::cross(vad, vab));
    double vd6 = glm::dot(vap, glm::cross(vab, vac));
    double v6 = 1 / glm::dot(vab, glm::cross(vac, vad));
    vector<double> b_coords = {va6*v6, vb6*v6, vc6*v6, vd6*v6};
    // Vertex v1 = vertices.at(cell.v_ids[0]);
    // Vertex v2 = vertices.at(cell.v_ids[1]);
    // Vertex v3 = vertices.at(cell.v_ids[2]);
    // Vertex v4 = vertices.at(cell.v_ids[3]);

    // glm::mat3 T;
    // T[0] = glm::vec3(v1.x - v4.x, v1.y - v4.y, v1.z - v4.z);
    // T[1] = glm::vec3(v2.x - v4.x, v2.y - v4.y, v2.z - v4.z);
    // T[2] = glm::vec3(v3.x - v4.x, v3.y - v4.y, v3.z - v4.z);
    
    // glm::vec3 rr4 = glm::vec3(p.x - v4.x, p.y - v4.y, p.z - v4.z);
    // glm::mat3 T_in = glm::inverse(T);
    // glm::vec3 barycentric_coords = T_in * rr4;
    // vector<double> b_coords = {barycentric_coords[0], barycentric_coords[1], barycentric_coords[2]};
    return b_coords;
}

void print_vector(string message, vector<double> v) {
    cout << message << ": ";
    for (int i = 0; i < v.size(); i++) {
        cout << v[i] << " ";
    }
    cout << endl;
}

void Mesh::setNeighboringCorners() {
    patches = extractPatches();
    for (int i = 0; i < patches.size(); i++) {
        setCornerNeighbors(patches[i]);
    }
    setStepSizes();
    setCornersUVWcoords();
}

void Mesh::setCornerNeighbors(vector<Face> patch) {
    for (int i = 0; i < patch.size(); i++) {
        for (int j = i+1; j < patch.size(); j++) {
            if (containsEdges(patch[i].v_ids, patch[j].v_ids)) {
                patch[i].neighbors.push_back(j);
                patch[j].neighbors.push_back(i);
            }
        }
    }
    vector<Edge> boundary_edges;
    for (int i = 0; i < patch.size(); i++) {
        for (int j = 0; j < patch[i].v_ids.size(); j++) {
            if (vertices.at(patch[i].v_ids[j]).isBoundaryVertex && vertices.at(patch[i].v_ids[(j+1)%patch[i].v_ids.size()]).isBoundaryVertex) {
                Edge e(patch[i].v_ids[j], patch[i].v_ids[(j+1)%patch[i].v_ids.size()], boundary_edges.size());
                if (!edgeInNeighbors(e, patch, i)) {
                    boundary_edges.push_back(e);
                }
            }
        }
    }
    for (int i = 0; i < boundary_edges.size(); i++) {
        for (int j = i + 1; j < boundary_edges.size(); j++) {
            if (boundary_edges[i].v1 == boundary_edges[j].v1 || boundary_edges[i].v2 == boundary_edges[j].v2 ||
                boundary_edges[i].v1 == boundary_edges[j].v2 || boundary_edges[i].v2 == boundary_edges[j].v1) {
                boundary_edges[i].neighbors.push_back(j);
                boundary_edges[j].neighbors.push_back(i);
            }
        }
    }
    int src_id = -1;
    vector<Edge> queue = boundary_edges;
    vector<int> v_ids;
    while (!queue.empty()) {
        Edge e = boundary_edges[queue.back().id];
        if (boundary_edges[e.id].isVisited) {
            queue.pop_back();
        } else {
            if (src_id == -1) {
                src_id = e.v1;
            }
            if (vertices.at(src_id).isCornerVertex) {
                v_ids.push_back(src_id);
            }
            boundary_edges[e.id].isVisited = true;
            bool neighborFound = false;
            for (int i = 0; i < e.neighbors.size(); i++) {
                if (!boundary_edges[e.neighbors[i]].isVisited) {
                    if (src_id == boundary_edges[e.neighbors[i]].v1) {
                        src_id = boundary_edges[e.neighbors[i]].v2;
                        queue.push_back(boundary_edges[e.neighbors[i]]);
                        neighborFound = true;
                    } else if (src_id == boundary_edges[e.neighbors[i]].v2) {
                        src_id = boundary_edges[e.neighbors[i]].v1;
                        queue.push_back(boundary_edges[e.neighbors[i]]);
                        neighborFound = true;
                    }
                }
            }
            if (!neighborFound && !v_ids.empty()) {
                src_id = -1;
                if (v_ids.size() >= 4) {
                    for (int i = 0; i < v_ids.size(); i++) {
                        Vertex& v1 = vertices.at(v_ids.at(i));
                        Vertex& v2 = vertices.at(v_ids.at((i+1)%v_ids.size()));
                        corner_edges.emplace_back(v1.id, v2.id, corner_edges.size());
                        v1.neighboring_corners.push_back(v2.id);
                        // v2.neighboring_corners.push_back(v1.id);
                    }
                }
                mesh_faces.emplace_back(v_ids);
                v_ids.clear();
            }
        }
    }
    for (int i = 0; i < corner_edges.size(); i++) {
        for (int j = 0; j < corner_edges.size(); j++) {
            if (i == j) {
                continue;
            }
            if ((corner_edges.at(i).v1 == corner_edges.at(j).v1 && corner_edges.at(i).v2 == corner_edges.at(j).v2)
            || (corner_edges.at(i).v1 == corner_edges.at(j).v2 && corner_edges.at(i).v2 == corner_edges.at(j).v1)) {
                corner_edges.erase(corner_edges.begin() + i);
                i = 0;
            }
        }
    }
}

void Mesh::setCornersUVWcoords() {
    vector<glm::vec3> directions = {
        glm::vec3(1, 0, 0),
        glm::vec3(0, 1, 0),
        glm::vec3(0, 0, 1),
    };
    for (int i = 0; i < corner_edges.size(); i++) {
        Vertex& v1 = vertices.at(corner_edges.at(i).v1);
        Vertex& v2 = vertices.at(corner_edges.at(i).v2);
        glm::vec3 v = glm::normalize(glm::vec3(v2.x - v1.x, v2.y - v1.y, v2.z - v1.z));
        // stepSizeX = stepSizeY = stepSizeZ = stepSize;
        if (round(abs(glm::dot(v, directions[0]))) == 1) {
            double p1 = floor(v1.x / stepSizeX) * stepSizeX;
            double p2 = ceil(v1.x / stepSizeX) * stepSizeX;
            double p3 = floor(v2.x / stepSizeX) * stepSizeX;
            double p4 = ceil(v2.x / stepSizeX) * stepSizeX;
            if (abs(p1 - v1.x) / abs(v2.x - v1.x) < 1 && abs(p1 - v2.x) / abs(v2.x - v1.x) < 1) {
                v1.u = p1;
            } else {
                v1.u = p2;
            }
            if (abs(p3 - v1.x) / abs(v2.x - v1.x) < 1 && abs(p3 - v2.x) / abs(v2.x - v1.x) < 1) {
                v2.u = p3;
            } else {
                v2.u = p4;
            }
        } else if (round(abs(glm::dot(v, directions[1]))) == 1) {
            double p1 = floor(v1.y / stepSizeY) * stepSizeY;
            double p2 = ceil(v1.y / stepSizeY) * stepSizeY;
            double p3 = floor(v2.y / stepSizeY) * stepSizeY;
            double p4 = ceil(v2.y / stepSizeY) * stepSizeY;
            if (abs(p1 - v1.y) / abs(v2.y - v1.y) < 1 && abs(p1 - v2.y) / abs(v2.y - v1.y) < 1) {
                v1.v = p1;
            } else {
                v1.v = p2;
            }
            if (abs(p3 - v1.y) / abs(v2.y - v1.y) < 1 && abs(p3 - v2.y) / abs(v2.y - v1.y) < 1) {
                v2.v = p3;
            } else {
                v2.v = p4;
            }
        } else if (round(abs(glm::dot(v, directions[2]))) == 1) {
            double p1 = floor(v1.z / stepSizeZ) * stepSizeZ;
            double p2 = ceil(v1.z / stepSizeZ) * stepSizeZ;
            double p3 = floor(v2.z / stepSizeZ) * stepSizeZ;
            double p4 = ceil(v2.z / stepSizeZ) * stepSizeZ;
            if (abs(p1 - v1.z) / abs(v2.z - v1.z) < 1 && abs(p1 - v2.z) / abs(v2.z - v1.z) < 1) {
                v1.w = p1;
            } else {
                v1.w = p2;
            }
            if (abs(p3 - v1.z) / abs(v2.z - v1.z) < 1 && abs(p3 - v2.z) / abs(v2.z - v1.z) < 1) {
                v2.w = p3;
            } else {
                v2.w = p4;
            }
        }
    }
    for (int i = 0; i < corner_vertices.size(); i++) {
        Vertex& v1 = vertices.at(corner_vertices.at(i));
        vector<double> p1 = {v1.u, v1.v, v1.w};
        if (isPointOutside(p1, v1)) {
            for (int j = 0; j < v1.neighboring_corners.size(); j++) {
                Vertex& v2 = vertices.at(v1.neighboring_corners.at(j));
                vector<double> p2 = {v2.u, v2.v, v2.w};
                if (isPointOutside(p2, v2)) {
                    continue;
                }
                glm::vec3 v = glm::normalize(glm::vec3(v2.x - v1.x, v2.y - v1.y, v2.z - v1.z));
                if (round(abs(glm::dot(v, directions[0]))) == 1) {
                    if (v1.v != v2.v) {
                        v1.v = v2.v;
                        alignNeighbors(v1, 1);
                    }
                    if (v1.w != v2.w) {
                        v1.w = v2.w;
                        alignNeighbors(v1, 2);
                    }
                }
                if (round(abs(glm::dot(v, directions[1]))) == 1) {
                    if (v1.u != v2.u) {
                        v1.u = v2.u;
                        alignNeighbors(v1, 0);
                    }
                    if (v1.w != v2.w) {
                        v1.w = v2.w;
                        alignNeighbors(v1, 2);
                    }
                }
                if (round(abs(glm::dot(v, directions[2]))) == 1) {
                    if (v1.u != v2.u) {
                        v1.u = v2.u;
                        alignNeighbors(v1, 0);
                    }
                    if (v1.v != v2.v) {
                        v1.v = v2.v;
                        alignNeighbors(v1, 1);
                    }
                }
            }
        }    
    }
    
    for (int i = 0; i < corner_edges.size(); i++) {
        Vertex& v1 = vertices.at(corner_edges.at(i).v1);
        Vertex& v2 = vertices.at(corner_edges.at(i).v2);

        glm::vec3 v = glm::normalize(glm::vec3(v2.x - v1.x, v2.y - v1.y, v2.z - v1.z));
        if (round(abs(glm::dot(v, directions[0]))) == 1) {
            if (v1.v != v2.v) {
                v1.v = v2.v;
                alignNeighbors(v1, 1);
                i = 0;
            }
            if (v1.w != v2.w) {
                v1.w = v2.w;
                alignNeighbors(v1, 2);
                i = 0;
            }
        } else if (round(abs(glm::dot(v, directions[1]))) == 1) {
            if (v1.u != v2.u) {
                v1.u = v2.u;
                alignNeighbors(v1, 0);
                i = 0;
            }
            if (v1.w != v2.w) {
                v1.w = v2.w;
                alignNeighbors(v1, 2);
                i = 0;
            }
        } else if (round(abs(glm::dot(v, directions[2]))) == 1) {
            if (v1.v != v2.v) {
                v1.v = v2.v;
                alignNeighbors(v1, 1);
                i = 0;
            }
            if (v1.u != v2.u) {
                v1.u = v2.u;
                alignNeighbors(v1, 0);
                i = 0;
            }
        }
    }
}

void Mesh::alignNeighbors(Vertex& v, int direction) {
    vector<glm::vec3> directions = {
        glm::vec3(1, 0, 0),
        glm::vec3(0, 1, 0),
        glm::vec3(0, 0, 1),
    };
    
    for (int i = 0; i < v.neighboring_corners.size(); i++) {
        Vertex& n = vertices.at(v.neighboring_corners.at(i));
        glm::vec3 d = glm::normalize(glm::vec3(v.x - n.x, v.y - n.y, v.z - n.z));
        if (round(abs(glm::dot(d, directions[direction]))) == 1) {
            continue;
        }
        if (direction == 0 && n.u != v.u) {
            n.u = v.u;
            alignNeighbors(n, direction);
        } else if (direction == 1 && n.v != v.v) {
            n.v = v.v;
            alignNeighbors(n, direction);
        } else if (direction == 2 && n.w != v.w) {
            n.w = v.w;
            alignNeighbors(n, direction);
        }
    }
}

bool Mesh::isPointOutside(vector<double> p, Vertex& v) {
    bool output = true;
    Vertex p1(p[0], p[1], p[2], -1);
    for (int i = 0; i < v.cells_ids.size(); i++) {
        vector<double> b_coords = getBarycentricCoordinates(p1, v.cells_ids.at(i));
        if (b_coords[0] >= 0 && b_coords[0] <= 1 && b_coords[1] >= 0 && b_coords[1] <= 1 && b_coords[2] >= 0 && b_coords[2] <= 1) {
            output = false;
            break;
        }
    }
    return output;
}

vector<Edge> Mesh::extractFacesFromPatch(vector<Face> patch) {
    
    for (int i = 0; i < patch.size(); i++) {
        for (int j = i+1; j < patch.size(); j++) {
            if (containsVertices(patch[i].v_ids, patch[j].v_ids)) {
                patch[i].neighbors.push_back(j);
                patch[j].neighbors.push_back(i);
            }
        }
    }
    vector<Edge> boundary_edges;
    for (int i = 0; i < patch.size(); i++) {
        if (!patch[i].isCornerFace) {
            continue;
        }
        for (int j = 0; j < patch[i].v_ids.size(); j++) {
            // if (vertices.at(patch[i].v_ids[j]).isCornerVertex && vertices.at(patch[i].v_ids[(j+1)%patch[i].v_ids.size()]).isBoundaryVertex) {
                Edge e(patch[i].v_ids[j], patch[i].v_ids[(j+1)%patch[i].v_ids.size()], boundary_edges.size());
                if (!edgeInNeighbors(e, patch, i)) {
                    vector<double> direction = normalizeVector(getDirectionVector(vertices.at(e.v1), vertices.at(e.v2)));
                    Vertex corner = vertices.at(e.v1);
                    if (vertices.at(e.v2).isCornerVertex) {
                        direction = normalizeVector(getDirectionVector(vertices.at(e.v2), vertices.at(e.v1)));
                        corner = vertices.at(e.v2);
                    }
                    int f_id = i;
                    Vertex initial = corner;
                    Edge traceEdge;
                    while (true) {
                        int id = vertices.size();
                        Vertex end(initial.x + (stepSize * direction[0]), initial.y + (stepSize * direction[1]), initial.z + (stepSize * direction[2]), id);
                        vertices.push_back(end);
                        boundary_edges.emplace_back(initial.id, end.id, boundary_edges.size());    
                        vector<double> barycentric = getBarycentricCoordinates(end, patch[f_id]);
                        if (barycentric[0] >= 0 && barycentric[1] >= 0 && barycentric[2] >= 0) {
                            initial = end;
                            continue;
                        } else {
                            bool foundFace = false;
                            for (int k = 0; k < patch[f_id].neighbors.size(); k++) {
                                barycentric = getBarycentricCoordinates(end, patch[patch[f_id].neighbors[k]]);
                                if (barycentric[0] >= 0 && barycentric[1] >= 0 && barycentric[2] >= 0) {
                                    initial = end;
                                    foundFace = true;
                                    f_id = patch[f_id].neighbors[k];
                                    break;
                                }
                            }
                            if (foundFace) {
                                continue;
                            }
                        }
                        // for (int l = 0; l < patch[f_id].v_ids.size(); l++) {
                        //     Edge boundary_edge(patch[i].v_ids[j], patch[i].v_ids[(j+1)%patch[i].v_ids.size()], -1);
                        //     if (!edgeInNeighbors(e, patch, l)) {
                        //         Vertex b_v1 = vertices.at(boundary_edge.v1);
                        //         Vertex b_v2 = vertices.at(boundary_edge.v2);
                        //     }
                        // }
                        break;
                    }
                    f_id = i;
                    initial = corner;
                    direction[0] = -direction[0];
                    direction[1] = -direction[1];
                    direction[2] = -direction[2];
                    while (true) {
                        int id = vertices.size();
                        Vertex end(initial.x + (stepSize * direction[0]), initial.y + (stepSize * direction[1]), initial.z + (stepSize * direction[2]), id);
                        vertices.push_back(end);
                        boundary_edges.emplace_back(initial.id, end.id, boundary_edges.size());    
                        vector<double> barycentric = getBarycentricCoordinates(end, patch[f_id]);
                        if (barycentric[0] >= 0 && barycentric[1] >= 0 && barycentric[2] >= 0) {
                            initial = end;
                            continue;
                        } else {
                            bool foundFace = false;
                            for (int k = 0; k < patch[f_id].neighbors.size(); k++) {
                                barycentric = getBarycentricCoordinates(end, patch[patch[f_id].neighbors[k]]);
                                if (barycentric[0] >= 0 && barycentric[1] >= 0 && barycentric[2] >= 0) {
                                    initial = end;
                                    foundFace = true;
                                    f_id = patch[f_id].neighbors[k];
                                    break;
                                }
                            }
                            if (foundFace) {
                                continue;
                            }
                        }
                        break;
                    }
                }
            // }
        }
    }
    return boundary_edges;
}

vector<double> Mesh::getBarycentricCoordinates(Vertex p, Face f) {
    Vertex a = vertices.at(f.v_ids[0]);
    Vertex b = vertices.at(f.v_ids[1]);
    Vertex c = vertices.at(f.v_ids[2]);

    vector<double> ab = getDirectionVector(a, b);
    vector<double> ac = getDirectionVector(a, c);
    vector<double> ca = getDirectionVector(c, a);
    vector<double> cp = getDirectionVector(c, p);
    vector<double> ap = getDirectionVector(a, p);
    vector<double> bc = getDirectionVector(b, c);
    vector<double> bp = getDirectionVector(b, p);

    vector<double> N = getCrossProduct(ab, ac);
    vector<double> u_vec = getCrossProduct(ca, cp);
    vector<double> v_vec = getCrossProduct(ab, ap);
    vector<double> w_vec = getCrossProduct(bc, bp);
    
    double n = getDotProduct(N, N);
    double u = getDotProduct(u_vec, N);
    double v = getDotProduct(v_vec, N);
    double w = getDotProduct(w_vec, N);
    vector<double> barycentric = {u, v, w, n};
    return barycentric;
}

vector<Face> Mesh::extractFaces() {
    vector<vector<Face>> patches = extractPatches();
    vector<Face> faces;
    for (int i = 0; i < patches.size(); i++) {
    }
    return faces;
}

bool Mesh::edgeInNeighbors(Edge e, vector<Face> patch, int id) {
    for (int i = 0; i < patch[id].neighbors.size(); i++) {
        if (count(patch[patch[id].neighbors[i]].v_ids.begin(), patch[patch[id].neighbors[i]].v_ids.end(), e.v1) && 
            count(patch[patch[id].neighbors[i]].v_ids.begin(), patch[patch[id].neighbors[i]].v_ids.end(), e.v2)) {
            return true;
        }
    }
    return false;  
}

bool Mesh::containsEdges(vector<int> v_ids1, vector<int> v_ids2) {
    for (int i = 0; i < v_ids1.size(); i++) {
        if (count(v_ids2.begin(), v_ids2.end(), v_ids1[i]) && count(v_ids2.begin(), v_ids2.end(), v_ids1[(i+1)%v_ids1.size()])) {
            return true;
        }
    }
    return false;
}

bool Mesh::containsVertices(vector<int> v_ids1, vector<int> v_ids2) {
    for (int i = 0; i < v_ids1.size(); i++) {
        if (count(v_ids2.begin(), v_ids2.end(), v_ids1[i])) {
            return true;
        }
    }
    return false;
}

void Mesh::computeNormal(int v_id0, int v_id1, int v_id2) {
    Vertex v0 = vertices.at(v_id0);
    Vertex v1 = vertices.at(v_id1);
    Vertex v2 = vertices.at(v_id2);

    vector<double> a = getDirectionVector(v0, v1);
    vector<double> b = getDirectionVector(v0, v2);

    vector<double> normal = normalizeVector(getCrossProduct(a, b));
    vertices.at(v_id0).normal_x = round(abs(normal[0])) || vertices.at(v_id0).normal_x;
    vertices.at(v_id0).normal_y = round(abs(normal[1])) || vertices.at(v_id0).normal_y;
    vertices.at(v_id0).normal_z = round(abs(normal[2])) || vertices.at(v_id0).normal_z;
    vertices.at(v_id0).normal_sum = vertices.at(v_id0).normal_x + vertices.at(v_id0).normal_y + vertices.at(v_id0).normal_z;
}

void Mesh::setStepSize(Face& f) {
    for (int j = 0; j < f.v_ids.size(); j++) {
        Vertex v1 = vertices.at(f.v_ids.at(j%f.v_ids.size()));
        Vertex v2 = vertices.at(f.v_ids.at((j+1)%f.v_ids.size()));
        double length = glm::distance(glm::vec3(v1.x, v1.y, v1.z), glm::vec3(v2.x, v2.y, v2.z));
        double factor = 1;
        if (length > 0 && stepSize > length * factor) {
            stepSize = length * factor;
        }
    }
}

void Mesh::setStepSizes() {
    vector<glm::vec3> directions = {
        glm::vec3(1, 0, 0),
        glm::vec3(0, 1, 0),
        glm::vec3(0, 0, 1),
    };
    for (int i = 0; i < corner_edges.size(); i++) {
        Vertex& v1 = vertices.at(corner_edges.at(i).v1);
        Vertex& v2 = vertices.at(corner_edges.at(i).v2);
        glm::vec3 v = glm::normalize(glm::vec3(v2.x - v1.x, v2.y - v1.y, v2.z - v1.z));
        double factor = 0.1;
        double length = glm::length(glm::vec3(v2.x - v1.x, v2.y - v1.y, v2.z - v1.z));
        if (round(abs(glm::dot(v, directions[0]))) == 1 && stepSizeX > (length * factor)) {
            stepSizeX = length * factor;
        } else if (round(abs(glm::dot(v, directions[1]))) == 1 && stepSizeY > (length * factor)) {
            stepSizeY = length * factor;
        } else if (round(abs(glm::dot(v, directions[2]))) == 1 && stepSizeZ > (length * factor)) {
            stepSizeZ = length * factor;
        }
    }
}

double Mesh::getVectorLength(vector<double> v) {
    return sqrt(pow((v[0]), 2) + pow((v[1]), 2) + pow((v[2]), 2));
}

vector<double> Mesh::getDirectionVector(Vertex v1, Vertex v2) {
    vector<double> d = {doubleDifference(v1.x, v2.x), doubleDifference(v1.y, v2.y), doubleDifference(v1.z, v2.z)};
    return d;
}

double Mesh::doubleDifference(double a, double b) {
    double diff = b - a;
    if (abs(diff) < 0.0001) {
        return 0;
    }
    return diff;
}

vector<double> Mesh::normalizeVector(vector<double> v) {
    double length = getVectorLength(v);
    vector<double> new_v = {v[0] / length, v[1] / length, v[2] / length};
    return new_v;
}

vector<double> Mesh::getCrossProduct(vector<double> a, vector<double> b) {
    vector<double> product = {(a[1] * b[2]) - (a[2] * b[1]), (a[2] * b[0]) - (a[0] * b[2]), (a[0] * b[1]) - (a[1] * b[0])};
    return product;
}

double Mesh::getDotProduct(vector<double> a, vector<double> b) {
    return (a[0] * b[0]) + (a[1] * b[1]) + (a[2] * b[2]);
}

vector<double> Mesh::rotateVector(vector<double> v, double angle, int axis) {
    angle = angle * PI / 180;
    glm::vec3 in = glm::vec3(v[0], v[1], v[2]);
    glm::mat3 rX;
    rX[0] = glm::vec3(1, 0, 0);
    rX[1] = glm::vec3(0, cos(angle), -sin(angle));
    rX[2] = glm::vec3(0, sin(angle), cos(angle));
    glm::mat3 rY;
    rY[0] = glm::vec3(cos(angle), 0, sin(angle));
    rY[1] = glm::vec3(0, 1, 0);
    rY[2] = glm::vec3(-sin(angle), 0, cos(angle));
    glm::mat3 rZ;
    rZ[0] = glm::vec3(cos(angle), -sin(angle), 0);
    rZ[1] = glm::vec3(sin(angle), cos(angle), 0);
    rZ[2] = glm::vec3(0, 0, 1);
    
    glm::vec3 r;
    if (axis == 0) {
        r = rX * in;
    } else if (axis == 1) {
        r = rY * in;
    } else if (axis == 2) {
        r = rZ * in;
    }
    r = glm::normalize(r);
    vector<double> rotation = {r[0], r[1], r[2]};
    return rotation;
}

void Mesh::scale() {
    double scalingFactor = 1 / stepSize;
    
    for (int i = 0; i < vertices.size(); i++) {
        Vertex& v = vertices.at(i);
        v.x = v.x * scalingFactor;
        v.y = v.y * scalingFactor;
        v.z = v.z * scalingFactor;
    }
}