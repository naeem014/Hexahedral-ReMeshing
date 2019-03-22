#include<iostream>
#include<math.h>
#include<algorithm>
#include<limits>
#include<glm/glm.hpp>

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

    // double a[] = {v1.x - v0.x, v1.y - v0.y, v1.z - v0.z};
    // double b[] = {v2.x - v0.x, v2.y - v0.y, v2.z - v0.z};

    vector<double> a = getDirectionVector(v0, v1);
    vector<double> b = getDirectionVector(v0, v2);

    vector<double> normal = normalizeVector(getCrossProduct(a, b));
    // normal[0] = fabs(normal[0]);
    // normal[1] = fabs(normal[1]);
    // normal[2] = fabs(normal[2]);
    // double normal_x = (a[1] * b[2]) - (a[2] * b[1]);
    // double normal_y = (a[2] * b[0]) - (a[0] * b[2]);
    // double normal_z = (a[0] * b[1]) - (a[1] * b[0]);

    // double length = sqrt((normal_x * normal_x) + (normal_y * normal_y) + (normal_z * normal_z));

    // vector<double> normal = {abs(normal_x / length), abs(normal_y / length), abs(normal_z / length)};
    return normal;
}

vector<Edge> Mesh::getBoundaryEdges() {
    vector<Edge> boundary_edges;
    for (int i = 0; i < corner_cells.size(); i++) {
        // Cell& c = corner_cells.at(3);
        Cell& c = corner_cells.at(i);
        for (int j = 0; j < c.faces.size(); j++) {
            if (c.faces.at(j).isCornerFace) {
                vector<Edge> edges;
                Face& f = c.faces.at(j);
                int corner_index = getCornerIndex(f);
                Vertex v0 = vertices.at(f.v_ids[corner_index]);
                Vertex v1 = vertices.at(f.v_ids[(corner_index + 1) % f.v_ids.size()]);
                Vertex v2 = vertices.at(f.v_ids[(corner_index + 2) % f.v_ids.size()]);
                vector<double> a = getDirectionVector(v0, v1);
                vector<double> b = getDirectionVector(v0, v2);
                // vector<double> direction = normalizeVector(getCrossProduct(a, b));
                vector<double> direction = getTraceDirection(c.faces.at(j));
                if (direction.empty()) {
                    continue;
                }
                // vector<double> direction = getPlaneNormal(c.faces.at(j));
                edges = traceLineFromCorner(v0, c.id, direction);
                boundary_edges.insert(boundary_edges.begin(), edges.begin(), edges.end());
                vector<double> directionN = {-direction[0], -direction[1], -direction[2]};
                edges = traceLineFromCorner(v0, c.id, directionN);
                boundary_edges.insert(boundary_edges.begin(), edges.begin(), edges.end());
            }
            // break; //remove
        }
    //     break; // remove
    }
    return boundary_edges;
}

vector<double> Mesh::getTraceDirection(Face& f) {
    vector<double> direction;
    int corner_index = getCornerIndex(f);
    
    Vertex v0 = vertices.at(f.v_ids[corner_index]);
    Vertex v1 = vertices.at(f.v_ids[(corner_index + 1) % f.v_ids.size()]);
    Vertex v2 = vertices.at(f.v_ids[(corner_index + 2) % f.v_ids.size()]);
    // cout << "v0: " << v0.x << " " << v0.y << " " << v0.z << endl;
    // cout << "v1: " << v1.x << " " << v1.y << " " << v1.z << endl;
    // cout << "v2: " << v2.x << " " << v2.y << " " << v2.z << endl;
    vector<double> a = normalizeVector(getDirectionVector(v0, v1));
    vector<double> b = normalizeVector(getDirectionVector(v0, v2));
    vector<double> face_normal = normalizeVector(getCrossProduct(a, b));
    // print_vector("a", a);
    // print_vector("b", b);
    for (int i = 0; i < 3; i++) {
        vector<double> n = normalizeVector(rotateVector(face_normal, 90, i));
        // print_vector("n",n);
        double dotA = fabs((getDotProduct(n, a)) / ((getVectorLength(n)) * (getVectorLength(a))));
        double dotB = fabs((getDotProduct(n, b)) / ((getVectorLength(n)) * (getVectorLength(b))));
        // cout << "dotA: " << dotA << endl;
        // cout << "dotB: " << dotB << endl;
        
        // cout << "==============================" << endl;
        if (dotA == 1) {
            // direction = normalizeVector(a);
            direction = a;
            break;
        } else if (dotB == 1) {
            // direction = normalizeVector(b);
            direction = b;
            break;
        }
    }
    return direction;
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

vector<Edge> Mesh::traceLineFromCorner(Vertex& corner, int c_id, vector<double> direction) {
    vector<Edge> boundary_edges;
    Cell& c = cells.at(c_id);
    Vertex initial = corner;
    while (true) {
        int id = vertices.size();
        Vertex end(initial.x + (stepSize * direction[0]), initial.y + (stepSize * direction[1]), initial.z + (stepSize * direction[2]), id);
        vertices.push_back(end);
        boundary_edges.emplace_back(initial.id, end.id, boundary_edges.size());    
        vector<double> b_coords = getBarycentricCoordinates(end, c_id);
        // cout << "a: " << b_coords[0] << " b: " << b_coords[1] << " c: " << b_coords[2] << " d: " << b_coords[3] << endl;
        // cout << "======================================================================================" << endl;
        if (b_coords[0] >= 0 && b_coords[1] >= 0 && b_coords[2] >= 0 && b_coords[3] >= 0) {
            initial = end;
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
                    if (b_coords[0] >= 0 && b_coords[1] >= 0 && b_coords[2] >= 0 && b_coords[3] >= 0) {
                        // cout << "a: " << b_coords[0] << " b: " << b_coords[1] << " c: " << b_coords[2] << " d: " << b_coords[3] << endl;
                        // cout << "||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||" << endl;
                    
                        initial = end;
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

    vector<double> ab = getDirectionVector(a, b);
    vector<double> ac = getDirectionVector(a, c);
    vector<double> ad = getDirectionVector(a, d);
    vector<double> ap = getDirectionVector(a, p);
    vector<double> bp = getDirectionVector(b, p);
    vector<double> bd = getDirectionVector(b, d);
    vector<double> bc = getDirectionVector(b, c);

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

    double V = 1 / fabs(getDotProduct(ad, getCrossProduct(ab, ac)));

    vector<double> b_coords = { getDotProduct(bp, getCrossProduct(bd, bc)) * V, 
                                getDotProduct(ap, getCrossProduct(ac, ad)) * V,
                                getDotProduct(ap, getCrossProduct(ad, ab)) * V,
                                getDotProduct(ap, getCrossProduct(ab, ac)) * V };

    return b_coords;
}

void print_vector(string message, vector<double> v) {
    cout << message << ": ";
    for (int i = 0; i < v.size(); i++) {
        cout << v[i] << " ";
    }
    cout << endl;
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
        vector<double> diff = getDirectionVector(v1, v2);
        double length = getVectorLength(diff);
        if (length > 0 && stepSize > length * 1) {
            stepSize = length * 1;
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
    vector<double> rotation = {v[0], v[1], v[2]};
    angle = angle * PI / 180;
    if (axis == 0) {
        rotation[1] = doubleDifference(v[2] * sin(angle), v[1] * cos(angle));
        rotation[2] = doubleDifference(v[2] * cos(angle), v[1] * sin(angle));
    } else if (axis == 1) {
        rotation[0] = doubleDifference(v[2] * sin(angle), v[0] * cos(angle));
        rotation[2] = doubleDifference(v[2] * cos(angle), -v[0] * sin(angle));
    } else if (axis == 2) {
        rotation[0] = doubleDifference(v[1] * sin(angle), v[0] * cos(angle));
        rotation[1] = doubleDifference(v[1] * cos(angle), v[0] * sin(angle));
    }
    rotation[0] = fabs(rotation[0]);
    rotation[1] = fabs(rotation[1]);
    rotation[2] = fabs(rotation[2]);
    return rotation;
}