#include<iostream>
#include<math.h>
#include<algorithm>
#include<limits>
#include "glm/glm.hpp"
#include "glm/gtx/rotate_vector.hpp"

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

    // vector<double> a = getDirectionVector(v0, v1);
    // vector<double> b = getDirectionVector(v0, v2);

    // vector<double> normal = normalizeVector(getCrossProduct(a, b));
    // normal[0] = fabs(normal[0]);
    // normal[1] = fabs(normal[1]);
    // normal[2] = fabs(normal[2]);
    double normal_x = (a[1] * b[2]) - (a[2] * b[1]);
    double normal_y = (a[2] * b[0]) - (a[0] * b[2]);
    double normal_z = (a[0] * b[1]) - (a[1] * b[0]);

    double length = sqrt((normal_x * normal_x) + (normal_y * normal_y) + (normal_z * normal_z));

    vector<double> normal = {abs(normal_x / length), abs(normal_y / length), abs(normal_z / length)};
    return normal;
}

vector<Edge> Mesh::getBoundaryEdges() {
    vector<Edge> boundary_edges;
    /*for (int i = 0; i < corner_vertices.size(); i++) {
        Vertex& c1 = vertices.at(corner_vertices.at(i));
        for (int j = 0; j < c1.neighboring_corners.size(); j++) {
            Vertex& c2 = vertices.at(c1.neighboring_corners.at(j));
            glm::vec3 V0 = glm::vec3(c1.x, c1.y, c1.z);
            glm::vec3 V1 = glm::vec3(c2.x, c2.y, c2.z);
            glm::vec3 d = glm::normalize(V1 - V0);
            vector<double> direction = {d[0], d[1], d[2]};
            for (int k = 0; k < corner_cells.size(); k++) {
                Cell& c = corner_cells.at(k);
                if (c.isCornerCell && count(c.v_ids.begin(), c.v_ids.end(), c1.id)) {
                    vector<Edge> edges = traceLineFromCorner(c1.id, c.id, direction);
            //         boundary_edges.insert(boundary_edges.begin(), edges.begin(), edges.end());
                }
            }
        }
    }*/
    for (int i = 3; i < corner_cells.size(); i++) {
        // Cell& c = corner_cells.at(3);
        Cell& c = corner_cells.at(i);
        for (int j = 0; j < c.v_ids.size(); j++) {
            if (!vertices.at(c.v_ids.at(j)).isCornerVertex) {
                continue;
            }
            // if (c.faces.at(j).isCornerFace) {
                vector<Edge> edges;
                // Face& f = c.faces.at(j);
                // int corner_index = getCornerIndex(f);
                Vertex& corner_v = vertices.at(c.v_ids[j]);
                // cout << "Before neighboring corners" << endl;
                // cout << corner_v.neighboring_corners.size() << endl;
                for (int k = 0; k < corner_v.neighboring_corners.size(); k++) {
                    // cout << "Inside neighboring corners loop: " << k << endl;
                    Vertex& corner_v2 = vertices.at(corner_v.neighboring_corners.at(k));
                    glm::vec3 V0 = glm::vec3(corner_v.x, corner_v.y, corner_v.z);
                    glm::vec3 V1 = glm::vec3(corner_v2.x, corner_v2.y, corner_v2.z);
                    // cout << V0[0] << " " << V0[1] << " " << V0[2] << endl;
                    // cout << V1[0] << " " << V1[1] << " " << V1[2] << endl;
                    // glm::vec3 d = V1 - V0;
                    glm::vec3 d = V1 - V0;
                    // glm::vec3 d = glm::normalize(V1 - V0);
                    vector<double> d1 = getDirectionVector(corner_v, corner_v2);
                    vector<double> direction1 = {d[0], d[1], d[2]};

                    // print_vector("direction1", direction1);
                    // cout << "Before tracing" << endl;
                    edges = traceLineFromCorner(corner_v.id, c.id, direction1);
                    // cout << "Corner: " << i << " vertex: " << j << " number of neighbors: " << corner_v.neighboring_corners.size() <<  " neighbor: " << k << endl;
                    boundary_edges.insert(boundary_edges.begin(), edges.begin(), edges.end());
                    // cout << "After concatenating edges" << endl;
                    break;
                }
                // cout << "Neighboring corners loop finished" << endl;
                // Vertex v0 = vertices.at(f.v_ids[corner_index]);
                // Vertex v1 = vertices.at(f.v_ids[(corner_index + 1) % f.v_ids.size()]);
                // Vertex v2 = vertices.at(f.v_ids[(corner_index + 2) % f.v_ids.size()]);
                // vector<double> a = getDirectionVector(v0, v1);
                // vector<double> b = getDirectionVector(v0, v2);
                // glm::vec3 V0 = glm::vec3(v0.x, v0.y, v0.z);
                // glm::vec3 V1 = glm::vec3(v1.x, v1.y, v1.z);
                // glm::vec3 V2 = glm::vec3(v2.x, v2.y, v2.z);
                // glm::vec3 a = V1 - V0;
                // glm::vec3 b = V2 - V0;
                // glm::vec3 calc = glm::normalize(glm::cross(a, b));
                // vector<double> direction1 = {calc[0], calc[1], calc[2]};
                // vector<double> direction = normalizeVector(getCrossProduct(a, b));
                // Edge e1(v0.id, v1.id, -1);
                // Edge e2(v0.id, v2.id, -1);
                // vector<double> direction1 = getTraceDirection(f, e1, c.id);
                // vector<double> direction2 = getTraceDirection(f, e2, c.id);
                // if (!direction1.empty()) {
                //     edges = traceLineFromCorner(v0, c.id, direction1);
                //     boundary_edges.insert(boundary_edges.begin(), edges.begin(), edges.end());
                    // vector<double> directionN1 = {-direction1[0], -direction1[1], -direction1[2]};
                    // edges = traceLineFromCorner(v0, c.id, directionN1);
                    // boundary_edges.insert(boundary_edges.begin(), edges.begin(), edges.end());
                // }
                // if (!direction2.empty()) {
                //     edges = traceLineFromCorner(v0, c.id, direction2);
                //     boundary_edges.insert(boundary_edges.begin(), edges.begin(), edges.end());
                //     vector<double> directionN2 = {-direction2[0], -direction2[1], -direction2[2]};
                //     edges = traceLineFromCorner(v0, c.id, directionN2);
                //     boundary_edges.insert(boundary_edges.begin(), edges.begin(), edges.end());
                // }
                // vector<double> direction = getPlaneNormal(c.faces.at(j));
            // }
            break; //remove
        }
        break; // remove
    }
    return boundary_edges;
}

vector<double> Mesh::getTraceDirection(Face& f, Edge e, int cell_id) {
    vector<double> direction;
    glm::vec3 d(vertices.at(e.v2).x - vertices.at(e.v1).x, vertices.at(e.v2).y - vertices.at(e.v1).y, vertices.at(e.v2).z - vertices.at(e.v1).z);
    Cell& c = cells.at(cell_id);
    bool foundFace = false;
    for (int i = 0; i < c.neighbors.size(); i++) {
        if (cells.at(c.neighbors.at(i)).isSurfaceCell) {
            Cell& neighbor_cell = cells.at(c.neighbors.at(i));
            for (int j = 0; j < neighbor_cell.faces.size(); j++) {
                Face& n_f = neighbor_cell.faces.at(j);
                if (n_f.isSurfaceFace && count(n_f.v_ids.begin(), n_f.v_ids.end(), e.v1) && count(n_f.v_ids.begin(), n_f.v_ids.end(), e.v2)) {
                    glm::vec3 a1(vertices.at(f.v_ids[0]).x - vertices.at(f.v_ids[1]).x, vertices.at(f.v_ids[0]).y - vertices.at(f.v_ids[1]).y, vertices.at(f.v_ids[0]).z - vertices.at(f.v_ids[1]).z);
                    glm::vec3 b1(vertices.at(f.v_ids[0]).x - vertices.at(f.v_ids[2]).x, vertices.at(f.v_ids[0]).y - vertices.at(f.v_ids[2]).y, vertices.at(f.v_ids[0]).z - vertices.at(f.v_ids[2]).z);
                    glm::vec3 a2(vertices.at(n_f.v_ids[0]).x - vertices.at(n_f.v_ids[1]).x, vertices.at(n_f.v_ids[0]).y - vertices.at(n_f.v_ids[1]).y, vertices.at(n_f.v_ids[0]).z - vertices.at(n_f.v_ids[1]).z);
                    glm::vec3 b2(vertices.at(n_f.v_ids[0]).x - vertices.at(n_f.v_ids[2]).x, vertices.at(n_f.v_ids[0]).y - vertices.at(n_f.v_ids[2]).y, vertices.at(n_f.v_ids[0]).z - vertices.at(n_f.v_ids[2]).z);
                    glm::vec3 c1 = glm::normalize(glm::cross(a1, b1));
                    glm::vec3 c2 = glm::normalize(glm::cross(a2, b2));
                    double dotP = abs(glm::dot(c1, c2) / (glm::length(c1) * glm::length(c2)));
                    if (dotP < 0.5) {
                        d = glm::normalize(d);
                        direction = {d[0], d[1], d[2]};
                        break;
                    }
                }
            }
            if (foundFace) {
                break;
            }
        }
    }
    return direction;
    
    // int corner_index = getCornerIndex(f);
    
    // Vertex v0 = vertices.at(f.v_ids[corner_index]);
    // Vertex v1 = vertices.at(f.v_ids[(corner_index + 1) % f.v_ids.size()]);
    // Vertex v2 = vertices.at(f.v_ids[(corner_index + 2) % f.v_ids.size()]);
    // cout << "v0: " << v0.x << " " << v0.y << " " << v0.z << endl;
    // cout << "v1: " << v1.x << " " << v1.y << " " << v1.z << endl;
    // cout << "v2: " << v2.x << " " << v2.y << " " << v2.z << endl;
    // glm::vec3 a = glm::vec3(v1.x, v1.y, v1.z) - glm::vec3(v0.x, v0.y, v0.z);
    // glm::vec3 b = glm::vec3(v2.x, v2.y, v2.z) - glm::vec3(v0.x, v0.y, v0.z);
    // glm::vec3 normal = glm::cross(a, b);
    // vector<double> a = normalizeVector(getDirectionVector(v0, v1));
    // vector<double> b = normalizeVector(getDirectionVector(v0, v2));
    // vector<double> face_normal = {normal[0], normal[1], normal[2]};
    // print_vector("a", a);
    // print_vector("b", b);
    // for (int i = 0; i < 3; i++) {
    //     vector<double> n_ = rotateVector(face_normal, 90, i);
    //     glm::vec3 n = glm::vec3(n_[0], n_[1], n_[2]);
    //     // print_vector("n",n);
    //     // double dotA = fabs((getDotProduct(n, a)) / ((getVectorLength(n)) * (getVectorLength(a))));
    //     // double dotB = fabs((getDotProduct(n, b)) / ((getVectorLength(n)) * (getVectorLength(b))));
    //     double dotA = glm::dot(n, a) / (glm::length(n) * glm::length(a));
    //     double dotB = glm::dot(n, b) / (glm::length(n) * glm::length(b));
    //     cout << "dotA: " << dotA << endl;
    //     cout << "dotB: " << dotB << endl;
        
    //     cout << "==============================" << endl;
    //     if (dotA == 1 || dotA == -1) {
    //         // direction = normalizeVector(a);
    //         a = glm::normalize(a);
    //         direction = {a[0], a[1], a[2]};
    //         break;
    //     } else if (dotB == 1 || dotB == -1) {
    //         // direction = normalizeVector(b);
    //         b = glm::normalize(b);
    //         direction = {b[0], b[1], b[2]};
    //         break;
    //     }
    // }
    // return direction;
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
    // Vertex initial = corner;
    stepSize = 0.001;
    glm::vec3 initial(corner.x, corner.y, corner.z);
    glm::vec3 d(stepSize * direction[0], stepSize * direction[1], stepSize * direction[2]);
    cout << d[0] << " " << d[1] << " " << d[2] << endl;
    // bool foundCell = false;
    while (true) {
        // cout << "Start of while loop" << endl;
        int id = vertices.size();
        // cout << "step size: " << stepSize << endl;
        // print_vector("direction: ", direction);
        // cout << "initial: " << initial.x << " " << initial.y << " " << initial.z << endl;
        // glm::vec3 new_end = glm::vec3(initial.x + (stepSize * direction[0]), initial.y + (stepSize * direction[1]), initial.z + (stepSize * direction[2]));
        Vertex start(initial[0], initial[1], initial[2], id);
        vertices.push_back(start);
        id = vertices.size();
        glm::vec3 new_end = initial + d;
        Vertex end(new_end[0], new_end[1], new_end[2], id);
        vertices.push_back(end);
        boundary_edges.emplace_back(start.id, end.id, boundary_edges.size());
        // if (foundCell) {
        //     stepSize = stepSize / 2;
        //     foundCell = false;
        // }
        // vertices.push_back(end);
        // boundary_edges.emplace_back(initial.id, end.id, boundary_edges.size());
        vector<double> b_coords = getBarycentricCoordinates(end, c_id);
        // cout << "C_ID: " << c_id << endl;
        cout << "a: " << b_coords[0] << " b: " << b_coords[1] << " c: " << b_coords[2] << endl;// << " d: " << b_coords[3] << endl;
        cout << b_coords[0] + b_coords[1] + b_coords[2] << endl;
        cout << c_id << endl;
        cout << "======================================================================================" << endl;
        if (b_coords[0] >= 0 && b_coords[0] <= 1 && b_coords[1] >= 0 && b_coords[1] <= 1 && b_coords[2] >= 0 && b_coords[2] <= 1) {// && b_coords[3] >= 0) {
            initial = new_end;
            continue;
        } else {
            // cout << "Entered neighbour condition" << endl;
            bool foundCell = false;
            // stepSize = stepSize * 2;
            for (int i = 0; i < c.v_ids.size(); i++) {
                // cout << "First for loop" << endl;
                Vertex& c_v = vertices.at(c.v_ids.at(i));
                for (int j = 0; j < c_v.cells_ids.size(); j++) {
                    if (c_v.cells_ids.at(j) == c_id) {
                        
                        // cout << "Continued " << j << " " << c_v.cells_ids.at(j) << endl;
                        continue;
                    }
                    // cout << "After if" << endl;
                    // int id_ = vertices.size();
                    // // cout << "id: " << id_ << endl;
                    // Vertex end_(initial.x + (2 * stepSize * direction[0]), initial.y + (2 * stepSize * direction[1]), initial.z + (2 * stepSize * direction[2]), -1);
                    // cout << "pushed edge" << endl;
                    // cout << end_.x << " " << end_.y << " " << end_.z << " " << end_.id  << endl;
                    // cout << j << " " << c_v.cells_ids.at(j) << endl;
                    b_coords = getBarycentricCoordinates(end, c_v.cells_ids.at(j));
                    // cout << "Got barycentric coordinates" << endl;
                    if (b_coords[0] >= 0 && b_coords[0] <= 1 && b_coords[1] >= 0 && b_coords[1] <= 1 && b_coords[2] >= 0 && b_coords[2] <= 1) {// && b_coords[3] >= 0) {
                        // cout << "a: " << b_coords[0] << " b: " << b_coords[1] << " c: " << b_coords[2] << endl;// << " d: " << b_coords[3] << endl;
                        // cout << "||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||" << endl;
                        // vertices.push_back(end_);
                        // boundary_edges.emplace_back(initial.id, end_.id, boundary_edges.size());    
                        // stepSize = stepSize / 2;
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
    // cout << "Inside bary func" << endl;
    Cell& cell = cells.at(cell_id);
    /*Vertex a = vertices.at(cell.v_ids[0]);
    Vertex b = vertices.at(cell.v_ids[1]);
    Vertex c = vertices.at(cell.v_ids[2]);
    Vertex d = vertices.at(cell.v_ids[3]);
    glm::vec3 ab = glm::vec3(b.x, b.y, b.z) - glm::vec3(a.x, a.y, a.z);
    glm::vec3 ac = glm::vec3(c.x, c.y, c.z) - glm::vec3(a.x, a.y, a.z);
    glm::vec3 ad = glm::vec3(d.x, d.y, d.z) - glm::vec3(a.x, a.y, a.z);
    glm::vec3 ap = glm::vec3(p.x, p.y, p.z) - glm::vec3(a.x, a.y, a.z);
    glm::vec3 bp = glm::vec3(p.x, p.y, p.z) - glm::vec3(b.x, b.y, b.z);
    glm::vec3 bd = glm::vec3(d.x, d.y, d.z) - glm::vec3(b.x, b.y, b.z);
    glm::vec3 bc = glm::vec3(c.x, c.y, c.z) - glm::vec3(b.x, b.y, b.z);
    // vector<double> ab = getDirectionVector(a, b);
    // vector<double> ac = getDirectionVector(a, c);
    // vector<double> ad = getDirectionVector(a, d);
    // vector<double> ap = getDirectionVector(a, p);
    // vector<double> bp = getDirectionVector(b, p);
    // vector<double> bd = getDirectionVector(b, d);
    // vector<double> bc = getDirectionVector(b, c);

    // cout << "a: " << a.x << " " << a.y << " " << a.z << endl;
    // cout << "b: " << b.x << " " << b.y << " " << b.z << endl;
    // cout << "c: " << c.x << " " << c.y << " " << c.z << endl;
    // cout << "d: " << d.x << " " << d.y << " " << d.z << endl;
    // cout << "p: " << p.x << " " << p.y << " " << p.z << endl;

    // cout << "----------------------------------------------------------------------" << endl;

    // print_vector("ab", ab);
    // print_vector("ac", ac);
    // print_vector("ad", ad);
    // print_vector("ap", ap);
    // print_vector("bp", bp);
    // print_vector("bd", bd);
    // print_vector("bc", bc);

    // // cout << "----------------------------------------------------------------------" << endl;

    // print_vector("bdXbc", getCrossProduct(bd, bc));
    // print_vector("acXad", getCrossProduct(ac, ad));
    // print_vector("adXab", getCrossProduct(ad, ab));
    // print_vector("abXac", getCrossProduct(ab, ac));

    // cout << "V: " << getDotProduct(ad, getCrossProduct(ab, ac));

    // cout << "----------------------------------------------------------------------" << endl;

    // double V = 1 / fabs(getDotProduct(ad, getCrossProduct(ab, ac)));
    double V = 1 / glm::dot(ab, glm::cross(ac, ad));
    // vector<double> b_coords = { getDotProduct(bp, getCrossProduct(bd, bc)) * V, 
    //                             getDotProduct(ap, getCrossProduct(ac, ad)) * V,
    //                             getDotProduct(ap, getCrossProduct(ad, ab)) * V,
    //                             getDotProduct(ap, getCrossProduct(ab, ac)) * V };
    double b1 = glm::dot(bp, glm::cross(bd, bc)) * V;
    double b2 = glm::dot(ap, glm::cross(ac, ad)) * V;
    double b3 = glm::dot(ap, glm::cross(ad, ab)) * V;
    double b4 = glm::dot(ap, glm::cross(ab, ac)) * V;
    // b1 = b1 > -0.00001 && b1 < 0.00001 ? 0 : b1;
    // b2 = b2 > -0.00001 && b2 < 0.00001 ? 0 : b2;
    // b3 = b3 > -0.00001 && b3 < 0.00001 ? 0 : b3;
    // b4 = b4 > -0.00001 && b4 < 0.00001? 0 : b4;
    vector<double> b_coords = { b1, b2, b3, b4};*/

    Vertex v1 = vertices.at(cell.v_ids[0]);
    Vertex v2 = vertices.at(cell.v_ids[1]);
    Vertex v3 = vertices.at(cell.v_ids[2]);
    Vertex v4 = vertices.at(cell.v_ids[3]);

    glm::mat3 T;
    T[0] = glm::vec3(v1.x - v4.x, v1.y - v4.y, v1.z - v4.z);
    T[1] = glm::vec3(v2.x - v4.x, v2.y - v4.y, v2.z - v4.z);
    T[2] = glm::vec3(v3.x - v4.x, v3.y - v4.y, v3.z - v4.z);
    glm::vec3 rr4 = glm::vec3(p.x - v4.x, p.y - v4.y, p.z - v4.z);
    glm::mat3 T_in = glm::inverse(T);
    glm::vec3 barycentric_coords = T_in * rr4;
    vector<double> b_coords = {barycentric_coords[0], barycentric_coords[1], barycentric_coords[2]};
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
    vector<vector<Face>> patches = extractPatches();
    // vector<int> c_vertices;
    for (int i = 0; i < patches.size(); i++) {
        setCornerNeighbors(patches[i]);
        // vector<int> cv = setCornerNeighbors(patches[i]);
        // c_vertices.insert(c_vertices.begin(), cv.begin(), cv.end());
        // break;
    }
    // return c_vertices;
}
/*
void Mesh::setCornerNeighbors(vector<Face> patch) {
    for (int i = 0; i < patch.size(); i++) {
        for (int j = i+1; j < patch.size(); j++) {
            if (containsVertices(patch[i].v_ids, patch[j].v_ids)) {
                patch[i].neighbors.push_back(j);
                patch[j].neighbors.push_back(i);
            }
        }
    }
    vector<Edge> boundary_edges;
    vector<Edge> queue;
    vector<int> c_vertices;
    for (int i = 0; i < patch.size(); i++) {
        if (!patch[i].isBoundaryFace) {
            continue;
        }
        for (int j = 0; j < patch[i].v_ids.size(); j++) {
            Edge e(patch[i].v_ids[j], patch[i].v_ids[(j+1)%patch[i].v_ids.size()], boundary_edges.size());
            if (!edgeInNeighbors(e, patch, i)) {
                boundary_edges.push_back(e);
                if (vertices.at(e.v1).isCornerVertex || vertices.at(e.v2).isCornerVertex) {
                    boundary_edges.back().isCornerEdge = true;
                    queue.push_back(boundary_edges.back());
                }
                break;
            }
        }
    }
    for (int i = 0; i < boundary_edges.size(); i++) {
        for (int j = i+1; j < boundary_edges.size(); j++) {
            if (boundary_edges[i].v1 == boundary_edges[j].v1 || boundary_edges[i].v1 == boundary_edges[j].v2
                || boundary_edges[i].v2 == boundary_edges[j].v1 || boundary_edges[i].v2 == boundary_edges[j].v2) {
                boundary_edges[i].neighbors.push_back(j);
                boundary_edges[j].neighbors.push_back(i);
            }
        }
    }
    while (!queue.empty()) {
        if (boundary_edges.at(queue.back().id).isVisited) {
            queue.pop_back();
        } else {
            boundary_edges.at(queue.back().id).isVisited = true;
            
            if (boundary_edges.at(queue.back().id).isCornerEdge) {
                Edge e = boundary_edges.at(queue.back().id);
                if (vertices.at(e.v1).isCornerVertex) {
                    c_vertices.push_back(e.v1);
                } else {
                    c_vertices.push_back(e.v2);
                }
            }
            for (int i = 0; i < boundary_edges.at(queue.back().id).neighbors.size(); i++) {
                if (!boundary_edges.at(boundary_edges.at(queue.back().id).neighbors.at(i)).isVisited) {
                    queue.push_back(boundary_edges.at(boundary_edges.at(queue.back().id).neighbors.at(i)));
                    break;
                }
            }
        }
    }
    // cout << c_vertices.size() << endl;
    for (int i = 0; i < c_vertices.size(); i++) {
        vertices.at(c_vertices[i]).neighboring_corners.push_back(vertices.at(c_vertices[(i+1)%c_vertices.size()]).id);
    }
    // return c_vertices;

        
        /*Face& f = patch[i];
        int corner_index = getCornerIndex(f);
        int corner_id = f.v_ids.at(corner_index);
        int v_id = corner_id;
        int patch_id = i;
        while (true) {
            for (int j = 0; j < patch[patch_id].v_ids.size(); j++) {
                Edge e(patch[patch_id].v_ids[j], patch[patch_id].v_ids[(j+1)%f.v_ids.size()], -1);
                if (!edgeInNeighbors(e, patch, patch_id)) {
                    if (v_id != e.v1) {
                        v_id = e.v1;
                    } else {
                        v_id = e.v2;
                    }
                    patch[patch_id].isVisited = true;
                    // cout << "corner id: " << corner_id << " " << "v id: " << v_id << endl;
                    break;
                }
            }
            if (vertices.at(v_id).isCornerVertex) {
                Vertex& v = vertices.at(v_id);
                Vertex& c = vertices.at(corner_id);
                if (v_id != corner_id && (v.neighbors.empty() || !count(v.neighbors.begin(), v.neighbors.end(), corner_id))) {
                    c.neighbors.push_back(v_id);
                }
                break;
            }
            for (int j = 0; j < patch[patch_id].neighbors.size(); j++) {
                // cout << "Checking neighbors" << endl;
                if (!patch[patch[patch_id].neighbors.at(j)].isVisited && count(patch[patch[patch_id].neighbors.at(j)].v_ids.begin(), patch[patch[patch_id].neighbors.at(j)].v_ids.end(), v_id)) {
                    patch_id = patch[patch_id].neighbors.at(j);
                    cout << patch_id << endl;
                    break;
                }
            }
        }
    }
}*/

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
    // vector<Face> faces;
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
                    // faces.emplace_back(v_ids);
                    for (int i = 0; i < v_ids.size(); i++) {
                        vertices.at(v_ids.at(i)).neighboring_corners.push_back(v_ids.at((i+1)%v_ids.size()));
                    }
                }
                v_ids.clear();
            }
        }
    }
}

vector<Edge> Mesh::extractFacesFromPatch(vector<Face> patch) {
    // setStepSize(patch);

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

/*vector<Face> Mesh::extractFacesFromPatch(vector<Face> patch) {
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
                if (!edgeInNeighbors(e, patch, patch[i].neighbors)) {
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
    vector<Face> faces;
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
                    faces.emplace_back(v_ids);
                }
                v_ids.clear();
            }
        }
    }
    return faces;
}*/

vector<Face> Mesh::extractFaces() {
    vector<vector<Face>> patches = extractPatches();
    vector<Face> faces;
    for (int i = 0; i < patches.size(); i++) {
        // vector<Face> f = extractFacesFromPatch(patches[i]);
        // faces.insert(faces.end(), f.begin(), f.end());
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
    // double normal_x = (a[1] * b[2]) - (a[2] * b[1]);
    // double normal_y = (a[2] * b[0]) - (a[0] * b[2]);
    // double normal_z = (a[0] * b[1]) - (a[1] * b[0]);

    // double length = sqrt((normal_x * normal_x) + (normal_y * normal_y) + (normal_z * normal_z));

    // double normal_array[] = {fabs(normal_x / length), fabs(normal_y / length), fabs(normal_z / length)};

    vertices.at(v_id0).normal_x = round(abs(normal[0])) || vertices.at(v_id0).normal_x;
    vertices.at(v_id0).normal_y = round(abs(normal[1])) || vertices.at(v_id0).normal_y;
    vertices.at(v_id0).normal_z = round(abs(normal[2])) || vertices.at(v_id0).normal_z;
    vertices.at(v_id0).normal_sum = vertices.at(v_id0).normal_x + vertices.at(v_id0).normal_y + vertices.at(v_id0).normal_z;
}

vector<Face> Mesh::refineFace(Face f) {
    vector<Face> faces;
    if (f.v_ids.size() == 4) {
        faces.push_back(f);
        return faces;
    }
    bool pointFound = false;
    for (int i = 0; i < f.v_ids.size(); i++) {
        Vertex p1 = vertices.at(f.v_ids.at(i));
        for (int j = 0; j < f.v_ids.size(); j++) {
            if (j == i || (j+1)%f.v_ids.size() == i || (j-1)%f.v_ids.size() == i || 
                (j+1)%f.v_ids.size() == (i-1)%f.v_ids.size() || (j-1)%f.v_ids.size() == (i+1)%f.v_ids.size()) {
                continue;
            }
            Vertex p2 = vertices.at(f.v_ids.at(j));
            Vertex p3 = vertices.at(f.v_ids.at((j+1)%f.v_ids.size()));
            vector<double> a = getDirectionVector(p2, p3);
            vector<double> b = getDirectionVector(p2, p1);
            double t = getDotProduct(a, b) / pow(getVectorLength(a), 2);
            if (t > 0 && t < 1) {
                cout << j << endl;
                cout << i << endl;
                cout << t << endl;
                vector<double> intersection = {p2.x + (t * a[0]), p2.y + (t * a[1]), p2.z + (t * a[2])};
                int id = vertices.size();
                vertices.emplace_back(intersection[0], intersection[1], intersection[2], id);
                pointFound = true;
                f.v_ids.insert(f.v_ids.begin() + j, id);
                break;    
            }
        }
        if (pointFound) {
            break;
        }
    }
    faces.push_back(f);
    return faces;
}

void Mesh::setStepSize(Face& f) {
    for (int j = 0; j < f.v_ids.size(); j++) {
        Vertex v1 = vertices.at(f.v_ids.at(j%f.v_ids.size()));
        Vertex v2 = vertices.at(f.v_ids.at((j+1)%f.v_ids.size()));
        double length = glm::distance(glm::vec3(v1.x, v1.y, v1.z), glm::vec3(v2.x, v2.y, v2.z));
        // vector<double> diff = getDirectionVector(v1, v2);
        // double length = getVectorLength(diff);
        if (length > 0 && stepSize > length * 0.5) {
            stepSize = length * 0.5;
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
    // vector<double> rotation = {v[0], v[1], v[2]};
    // if (axis == 0) {
    //     rotation[1] = doubleDifference(v[2] * sin(angle), v[1] * cos(angle));
    //     rotation[2] = doubleDifference(v[2] * cos(angle), v[1] * sin(angle));
    // } else if (axis == 1) {
    //     rotation[0] = doubleDifference(v[2] * sin(angle), v[0] * cos(angle));
    //     rotation[2] = doubleDifference(v[2] * cos(angle), -v[0] * sin(angle));
    // } else if (axis == 2) {
    //     rotation[0] = doubleDifference(v[1] * sin(angle), v[0] * cos(angle));
    //     rotation[1] = doubleDifference(v[1] * cos(angle), v[0] * sin(angle));
    // }
    // rotation[0] = fabs(rotation[0]);
    // rotation[1] = fabs(rotation[1]);
    // rotation[2] = fabs(rotation[2]);
    r = glm::normalize(r);
    vector<double> rotation = {r[0], r[1], r[2]};
    return rotation;
}