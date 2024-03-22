#include "halfEdgeMesh.h"
#include <iostream>
#include <algorithm> // For std::find_if

HalfEdgeMesh::HalfEdgeMesh(const std::vector<Eigen::Vector3f>& _vertices,
                           const std::vector<Eigen::Vector3i>& _faces) {
    buildHalfEdgeStructure(_vertices, _faces);
}

HalfEdgeMesh::~HalfEdgeMesh() {
    // Cleanup if needed, but usually, lists handle their own cleanup
}

void HalfEdgeMesh::buildHalfEdgeStructure(const std::vector<Eigen::Vector3f>& _vertices,
                                          const std::vector<Eigen::Vector3i>& _faces) {
    // Create vertices
    for (const auto& pos : _vertices) {
        vertices.emplace_back(pos);
    }

    // Create a map for edge lookup
    std::unordered_map<std::pair<int, int>, HalfEdge*, pair_hash> edgeMap;

    // Temporary storage to keep track of newly created halfEdges for easy access
    std::vector<HalfEdge*> tempHalfEdges;

    size_t faceIndex = 0;
    for (const auto& face : _faces) {
        HalfEdge* firstHalfEdge = nullptr;
        HalfEdge* prevHalfEdge = nullptr;

        for (int i = 0; i < 3; ++i) {
            int vStart = face[i];
            int vEnd = face[(i + 1) % 3];

            // Create halfEdge
            halfEdges.emplace_back();
            HalfEdge* newHalfEdge = &halfEdges.back();
            tempHalfEdges.push_back(newHalfEdge);

            // Set vertex
            auto vertexIt = std::next(vertices.begin(), vStart);
            newHalfEdge->vertex = &(*vertexIt);

            // Link halfEdges
            if (firstHalfEdge == nullptr) {
                firstHalfEdge = newHalfEdge;
            }
            if (prevHalfEdge != nullptr) {
                prevHalfEdge->next = newHalfEdge;
            }
            prevHalfEdge = newHalfEdge;

            // Set the vertex's halfEdge to one of its outgoing halfEdges
            if (vertexIt->halfEdge == nullptr) {
                vertexIt->halfEdge = newHalfEdge;
            }

            // Edge map handling
            std::pair<int, int> edgeKey(vStart, vEnd);
            edgeMap[edgeKey] = newHalfEdge;
        }

        // Close the loop
        prevHalfEdge->next = firstHalfEdge;

        // Setup twins
        for (int i = 0; i < 3; ++i) {
            int vStart = face[i];
            int vEnd = face[(i + 1) % 3];

            std::pair<int, int> edgeKey(vEnd, vStart);
            if (edgeMap.find(edgeKey) != edgeMap.end()) {
                HalfEdge* twinHalfEdge = edgeMap[edgeKey];
                tempHalfEdges[faceIndex * 3 + i]->twin = twinHalfEdge;
                twinHalfEdge->twin = tempHalfEdges[faceIndex * 3 + i];
            }
        }

        ++faceIndex;
    }
}

void HalfEdgeMesh::convertToMeshFormat(std::vector<Eigen::Vector3f>& outVertices,
                                       std::vector<Eigen::Vector3i>& outFaces) {
    // Clear output vectors
    outVertices.clear();
    outFaces.clear();

    // Fill outVertices
    for (const auto& vertex : vertices) {
        outVertices.push_back(vertex.position);
    }

    // An unordered set to keep track of visited halfEdges
    std::unordered_set<HalfEdge*> visited;

    for (auto& he : halfEdges) {
        if (visited.find(&he) == visited.end()) {
            std::vector<int> faceIndices;
            HalfEdge* start = &he;
            do {
                auto it = std::find_if(vertices.begin(), vertices.end(),
                                       [&start](const Vertex& v) { return &v == start->vertex; });
                int vertexIndex = std::distance(vertices.begin(), it);
                faceIndices.push_back(vertexIndex);

                visited.insert(start);
                start = start->next;
            } while (start != &he);

            if (faceIndices.size() == 3) {
                outFaces.push_back(Eigen::Vector3i(faceIndices[0], faceIndices[1], faceIndices[2]));
            }
        }
    }
}

//*************************************** Tests ***************************************//
#include <cassert>
// The version use assert
//bool HalfEdgeMesh::validate() const {
//    //***** Check for existence *****//
//    for (const Vertex& v : vertices) {
//        assert(v.halfEdge != nullptr && "Vertex without a halfEdge found.");
//    }

//    for (const HalfEdge& he : halfEdges) {
//        assert(he.next != nullptr && "An edge without a next found after collapsing.");
//    }

//    for (const HalfEdge& he : halfEdges) {
//        assert(he.twin != nullptr && "An edge without a twin found after collapsing.");
//    }

//    for (const HalfEdge& he : halfEdges) {
//        assert(he.vertex != nullptr && "An edge without a vertex found after collapsing.");
//    }

//    //***** Check for consistency *****//
//    for (const HalfEdge& he : halfEdges) {
//        assert(he.twin != nullptr && he.twin->twin == &he && "Twin HalfEdge inconsistency found.");
//    }

//    for (const HalfEdge& he : halfEdges) {
//        const HalfEdge* start = &he;
//        const HalfEdge* next = he.next;
//        while (next != nullptr && next != start) {
//            next = next->next;
//        }
//        assert(next == start && "Next HalfEdge does not form a loop.");
//    }

//    // Validate Face Consistency
//    for (const HalfEdge& he : halfEdges) {
//        const HalfEdge* start = &he;
//        const HalfEdge* current = he.next;
//        int count = 1;
//        while (current != nullptr && current != start && count <= 3) {
//            current = current->next;
//            count++;
//        }
//        assert(count >= 3 && "Face with less than 3 edges found.");
//    }

//    //***** Check for linkage *****//
//    for (const Vertex& v : vertices) {
//        const HalfEdge* start = v.halfEdge;
//        const HalfEdge* current = start;
//        do {
//            current = current->next ? current->next->twin : nullptr;
//        } while (current && current != start);

//        assert(current == start && "Vertex-HalfEdge linkage broken.");
//    }

//    // Validate HalfEdge Vertex Reference:
//    for (const HalfEdge& he : halfEdges) {
//        auto it = std::find_if(vertices.begin(), vertices.end(),
//                               [&he](const Vertex& v) { return &v == he.vertex; });
//        assert(it != vertices.end() && "HalfEdge points to a vertex not in the mesh.");
//    }

//    // Check for Duplicate HalfEdges:
//    for (auto it1 = halfEdges.begin(); it1 != halfEdges.end(); ++it1) {
//        for (auto it2 = std::next(it1); it2 != halfEdges.end(); ++it2) {
//            assert(!(it1->vertex == it2->vertex && it1->next->vertex == it2->next->vertex) && "Duplicate half-edges found.");
//        }
//    }

//    // Check for Isolated Half-Edges
//    for (const HalfEdge& he : halfEdges) {
//        assert(he.next != nullptr && he.twin != nullptr && "Isolated half-edge found.");
//    }

//    // Manifoldness Verification
//    std::unordered_set<std::pair<const Vertex*, const Vertex*>, pair_hash> edgeSet;
//    for (const HalfEdge& he : halfEdges) {
//        auto edge = std::make_pair(he.vertex, he.next->vertex);
//        assert(edgeSet.insert(edge).second && "Non-manifold edge detected.");
//    }

//    // Note: Check for Vertex Degree and Isolated Vertices, and Face Orientation and Consistency have been integrated earlier.

//    // If all tests pass
//    std::cout << "All tests passed" << std::endl;
//    return true; // This return statement could be removed since asserts handle the validation.
//}

// The version not use assert
bool HalfEdgeMesh::validate() const {
//    std::cout << "Start validating" << std::endl;
    //***** Check for exsistence *****//
    // Verify every vertex have a HalfEdge
    for (const Vertex& v : vertices) {
        if (v.halfEdge == nullptr) {
            std::cerr << "Vertex without a halfEdge found.\n";
            return false;
        }
    }

    // Verify every edge has a 'next'
    for (const HalfEdge& he : halfEdges) {
        if (!he.next) {
            std::cerr << "An edge without a next found after collapsing.\n";
            return false; // Indicates the function failed due to incorrect mesh state
        }
    }

    // Verify every edge has a 'twin'
    for (const HalfEdge& he : halfEdges) {
        if (!he.twin) {
            std::cerr << "An edge without a twin found after collapsing.\n";
            return false; // Indicates the function failed due to incorrect mesh state
        }
    }

    // Verify every edge has a 'vertex'
    for (const HalfEdge& he : halfEdges) {
        if (!he.vertex) {
            std::cerr << "An edge without a vertex found after collapsing.\n";
            return false; // Indicates the function failed due to incorrect mesh state
        }
    }

    //***** Check for consistency *****//
    // Twin HalfEdge Consistency:
    // Verify each half-edge has a twin and that the twin's twin is the original half-edge.
    for (const HalfEdge& he : halfEdges) {
        if (he.twin == nullptr || he.twin->twin != &he) {
            std::cerr << "Twin HalfEdge inconsistency found.\n";
            return false;
        }
    }

    // Next HalfEdge Consistency:
    // Ensure the next pointer of each half-edge eventually loops back to form a valid face.
    for (const HalfEdge& he : halfEdges) {
        const HalfEdge* start = &he;
        const HalfEdge* next = he.next;
        while (next != nullptr && next != start) {
            next = next->next;
        }
        if (next != start) {
            std::cerr << "Next HalfEdge does not form a loop.\n";
            return false;
        }
    }

    // Validate Face Consistency
    for (const HalfEdge& he : halfEdges) {
        const HalfEdge* start = &he;
        const HalfEdge* current = he.next;
        int count = 1; // Start has already been counted.
        while (current != nullptr && current != start && count <= 3) {
            current = current->next;
            count++;
        }
        if (count < 3) {
            std::cerr << "Face with less than 3 edges found.\n";
            return false;
        }
    }

    //***** Check for linkage *****//
    // Ensure Correct Vertex-HalfEdge Linkage
    for (const Vertex& v : vertices) {
        const HalfEdge* start = v.halfEdge;
        const HalfEdge* current = start;
        do {
            current = current->next->twin;
        } while (current && current != start);

        if (current != start) {
            std::cerr << "Vertex-HalfEdge linkage broken.\n";
            return false;
        }
    }

    // Validate HalfEdge Vertex Reference:
    // Ensure each half-edge points to a vertex present in the mesh.
    for (const HalfEdge& he : halfEdges) {
        if (std::find_if(vertices.begin(), vertices.end(),
                         [&he](const Vertex& v) { return &v == he.vertex; }) == vertices.end()) {
            std::cerr << "HalfEdge points to a vertex not in the mesh.\n";
            return false;
        }
    }

    // Check for Duplicate HalfEdges:
    // Ensure there are no duplicate half-edges (with the same vertex start and end points).
    for (auto it1 = halfEdges.begin(); it1 != halfEdges.end(); ++it1) {
        for (auto it2 = std::next(it1); it2 != halfEdges.end(); ++it2) {
            if (it1->vertex == it2->vertex && it1->next->vertex == it2->next->vertex) {
                std::cerr << "Duplicate half-edges found.\n";
                return false;
            }
        }
    }

    // Check for Isolated Half-Edges
    for (const HalfEdge& he : halfEdges) {
        if (he.next == nullptr || he.twin == nullptr) {
            std::cerr << "Isolated half-edge found.\n";
            return false;
        }
    }

    // Manifoldness Verification
    std::unordered_set<std::pair<const Vertex*, const Vertex*>, pair_hash> edgeSet;
    for (const HalfEdge& he : halfEdges) {
        auto edge = std::make_pair(he.vertex, he.next->vertex);
        if (!edgeSet.insert(edge).second) {
            std::cerr << "Non-manifold edge detected.\n";
            return false;
        }
    }

    // Check for Vertex Degree and Isolated Vertices
    for (const Vertex& v : vertices) {
        const HalfEdge* start = v.halfEdge;
        const HalfEdge* current = start;
        bool foundConnectedFace = false;

        do {
            if (!current) break; // Break if we somehow reach a null half-edge.
            if (current->next && current->next->vertex != &v) {
                foundConnectedFace = true;
                break; // Found a face connected to this vertex.
            }
            current = current->next ? current->next->twin : nullptr;
        } while (current != start);

        if (!foundConnectedFace) {
            std::cerr << "Isolated vertex or vertex without adequate degree found.\n";
            return false;
        }
    }

    // Check for Face Orientation and Consistency
    for (const HalfEdge& he : halfEdges) {
        const HalfEdge* start = &he;
        const HalfEdge* current = he.next;
        int loopCounter = 0;

        while (current && current != start) {
            current = current->next;
            loopCounter++;
            if (loopCounter > halfEdges.size()) { // Safety check to prevent infinite loops.
                std::cerr << "Face loop does not close properly.\n";
                return false;
            }
        }

        if (!current) {
            std::cerr << "Open face loop detected, indicating inconsistent face orientation.\n";
            return false;
        }
    }

    // Next HalfEdge Consistency:
    // Ensure the next pointer of each half-edge does not point to itself and eventually loops back.
    for (const HalfEdge& he : halfEdges) {
        const HalfEdge* current = &he;
        std::unordered_set<const HalfEdge*> visited;
        int stepCount = 0;

        do {
            if (visited.find(current) != visited.end()) {
                std::cerr << "Infinite loop detected in half-edge next pointers.\n";
                return false;
            }
            visited.insert(current);

            if (current->next == current) {
                std::cerr << "HalfEdge next pointer points to itself.\n";
                return false;
            }

            current = current->next;
            stepCount++;
            if (stepCount > halfEdges.size()) {
                std::cerr << "Excessive steps taken in half-edge next traversal, potential infinite loop.\n";
                return false;
            }
        } while (current != &he);
    }

    // Ensure that from every vertex, using twin->next can loop back to the original half-edge
    for (const Vertex& vertex : vertices) {
        const HalfEdge* startEdge = vertex.halfEdge;
        const HalfEdge* currentEdge = startEdge;
        std::unordered_set<const HalfEdge*> visitedEdges;

        do {
            // Check for the twin and next to avoid null pointer dereference
            if (!currentEdge->twin || !currentEdge->twin->next) {
                std::cerr << "Broken twin->next linkage.\n";
                return false;
            }

            // Move to the next half-edge around the vertex using twin->next
            currentEdge = currentEdge->twin->next;

            // Check if we have already visited this half-edge to detect a loop
            if (visitedEdges.find(currentEdge) != visitedEdges.end()) {
                std::cerr << "Infinite loop detected around a vertex.\n";
                return false;
            }
            visitedEdges.insert(currentEdge);

            // Also, ensure that we do not exceed the total number of half-edges in the mesh
            if (visitedEdges.size() > halfEdges.size()) {
                std::cerr << "Excessive half-edges visited around a vertex, potential infinite loop.\n";
                return false;
            }

        } while (currentEdge != startEdge);

        // If we have looped back to the starting edge, the linkage around this vertex is correct
    }

    // If all tests pass
    std::cout << "All tests passed" << std::endl;
    return true;
}

