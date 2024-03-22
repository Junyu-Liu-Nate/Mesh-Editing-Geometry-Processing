#include "operator.h"
#include <iostream>
#include <set>
#include <Eigen/Core>
#include <Eigen/Geometry>
#include <Eigen/SVD>
#include <random>
#include <queue>

Operator::Operator()
{

}

/********** 3 local operators **********/
bool Operator::flipEdge(HalfEdgeMesh& heMesh, HalfEdge* edge) {
    // Edge cannot be flipped if it does not have a twin or if it's a boundary edge.
    if (edge == nullptr || edge->twin == nullptr || !checkFlip(edge)) {
        return false;
    }

    // Define the half-edges involved in the operation
    HalfEdge* h0 = edge;
    HalfEdge* h1 = edge->next;
    HalfEdge* h2 = h1->next;
    HalfEdge* h3 = edge->twin;
    HalfEdge* h4 = h3->next;
    HalfEdge* h5 = h4->next;

    // Define the vertices involved
    Vertex* v0 = h0->vertex;
    Vertex* v1 = h2->vertex;
    Vertex* v2 = h1->vertex;
    Vertex* v3 = h5->vertex;

    // Edges h0 and h3 are the ones being flipped.

    // Re-assign half-edges' vertex
    h0->vertex = v1;
    h3->vertex = v3;

    // Re-assign vertices' half-edge
    v0->halfEdge = h4;
    v1->halfEdge = h0;
    v2->halfEdge = h1;
    v3->halfEdge = h3;

    // Re-assign half-edges' next pointers
    h0->next = h5;
    h5->next = h1;
    h1->next = h0;

    h3->next = h2;
    h2->next = h4;
    h4->next = h3;

    // Success
    return true;
}

std::pair<HalfEdge*, HalfEdge*> Operator::flipEdgeReturn(HalfEdgeMesh& heMesh, HalfEdge* edge) {
    // Edge cannot be flipped if it does not have a twin or if it's a boundary edge.
    if (edge == nullptr || edge->twin == nullptr || !checkFlip(edge)) {
        return {nullptr, nullptr}; // Adjusted to return a pair of nullptrs for failure cases
    }

    // Define the half-edges involved in the operation
    HalfEdge* h0 = edge;
    HalfEdge* h1 = edge->next;
    HalfEdge* h2 = h1->next;
    HalfEdge* h3 = edge->twin;
    HalfEdge* h4 = h3->next;
    HalfEdge* h5 = h4->next;

    // Define the vertices involved
    Vertex* v0 = h0->vertex;
    Vertex* v1 = h2->vertex;
    Vertex* v2 = h1->vertex;
    Vertex* v3 = h5->vertex;

    // Edges h0 and h3 are the ones being flipped.

    // Re-assign half-edges' vertex
    h0->vertex = v1;
    h3->vertex = v3;

    // Re-assign vertices' half-edge
    if (v0->halfEdge == h0) v0->halfEdge = h4; // Ensure only change if it points to the edge being flipped
    if (v1->halfEdge == h3) v1->halfEdge = h0;
    v2->halfEdge = h1;
    v3->halfEdge = h3;

    // Re-assign half-edges' next pointers
    h0->next = h5;
    h5->next = h1;
    h1->next = h0;

    h3->next = h2;
    h2->next = h4;
    h4->next = h3;

    // Return the modified edges
    return {h0, h3};
}

bool Operator::checkFlip(HalfEdge* edge) {
    Vertex* v1 = edge->vertex;
    Vertex* v2 = edge->twin->vertex;

    if (vertexDegree(v1) == 3 || vertexDegree(v2) == 3) {
        return false;
    }
    else {
        return true;
    }
}

HalfEdge* Operator::findFlipableEdge(HalfEdgeMesh& mesh) {
    for (auto& edge : mesh.halfEdges) {
        if (edge.twin != nullptr) { // Ensure it's not a boundary edge
            // Check if the edge forms part of a quadrilateral
            HalfEdge* next = edge.next;
            HalfEdge* twinNext = edge.twin->next;

            if (next != nullptr && twinNext != nullptr &&
                next->next != nullptr && twinNext->next != nullptr &&
                next->next->next == &edge && // Check if the edge loop is closed properly
                twinNext->next->next == edge.twin) { // Check twin's loop closure
                return &edge;
            }
        }
    }
    return nullptr; // No suitable edge found
}

std::vector<HalfEdge*> Operator::splitEdge(HalfEdgeMesh& heMesh, HalfEdge* edge) {
//    if (edge == nullptr || edge->twin == nullptr) {
//        return false;
//    }

    // Define the half-edges involved in the operation
    HalfEdge* h0 = edge;
    HalfEdge* h1 = edge->next;
    HalfEdge* h2 = h1->next;
    HalfEdge* h3 = edge->twin;
    HalfEdge* h4 = h3->next;
    HalfEdge* h5 = h4->next;

    // Define the vertices involved
    Vertex* v0 = h0->vertex;
    Vertex* v1 = h3->vertex; // This should be the vertex at the start of the twin edge
    Vertex* v2 = h2->vertex;
    Vertex* v3 = h5->vertex;

    // Calculate midpoint
    Eigen::Vector3f midpoint = (v0->position + v1->position) * 0.5f;
//    Eigen::Vector3f midpoint = 0.375*v0->position + 0.375*v1->position + 0.125*v2->position + 0.125*v3->position;

    // Create new vertex
    heMesh.vertices.emplace_back(midpoint);
    Vertex* newVertexPtr = &heMesh.vertices.back();

    // Create new half-edges and add them to the mesh's half-edges
    auto newHE1It = heMesh.halfEdges.emplace(heMesh.halfEdges.end());
    auto newHE2It = heMesh.halfEdges.emplace(heMesh.halfEdges.end());
    auto newHE3It = heMesh.halfEdges.emplace(heMesh.halfEdges.end());
    auto newHE4It = heMesh.halfEdges.emplace(heMesh.halfEdges.end());
    auto newHE5It = heMesh.halfEdges.emplace(heMesh.halfEdges.end());
    auto newHE6It = heMesh.halfEdges.emplace(heMesh.halfEdges.end());

    // Dereference iterators to get pointers to the actual half-edge elements
    HalfEdge* newHE1Ptr = &(*newHE1It);
    HalfEdge* newHE2Ptr = &(*newHE2It);
    HalfEdge* newHE3Ptr = &(*newHE3It);
    HalfEdge* newHE4Ptr = &(*newHE4It);
    HalfEdge* newHE5Ptr = &(*newHE5It);
    HalfEdge* newHE6Ptr = &(*newHE6It);

    // Update half-edges' vertices
    newHE1Ptr->vertex = newVertexPtr;
    newHE2Ptr->vertex = v2;
    newHE3Ptr->vertex = newVertexPtr;
    newHE4Ptr->vertex = newVertexPtr;
    newHE5Ptr->vertex = v3;
    newHE6Ptr->vertex = newVertexPtr;

    // Update half-edges' twins
    h0->twin = newHE4Ptr;
    newHE4Ptr->twin = h0;

    h3->twin = newHE1Ptr;
    newHE1Ptr->twin = h3;

    newHE2Ptr->twin = newHE3Ptr;
    newHE3Ptr->twin = newHE2Ptr;

    newHE6Ptr->twin = newHE5Ptr;
    newHE5Ptr->twin = newHE6Ptr;

    // Update original and new half-edges' next pointers
    h0->next = newHE3Ptr;
    newHE3Ptr->next = h2;
    h2->next = h0;

    h1->next = newHE2Ptr;
    newHE2Ptr->next = newHE1Ptr;
    newHE1Ptr->next = h1;

    h3->next = newHE6Ptr;
    newHE6Ptr->next = h5;
    h5->next = h3;

    h4->next = newHE5Ptr;
    newHE5Ptr->next = newHE4Ptr;
    newHE4Ptr->next = h4;

    // Update vertices' half-edge pointers
    newVertexPtr->halfEdge = newHE1Ptr;
    v0->halfEdge = h0;
    v1->halfEdge = h1;
    v2->halfEdge = h2;
    v3->halfEdge = h5;

    std::vector<HalfEdge*> newHalfEdges;
    newHalfEdges.push_back(newHE2Ptr);
    newHalfEdges.push_back(newHE3Ptr);
    newHalfEdges.push_back(newHE5Ptr);
    newHalfEdges.push_back(newHE6Ptr);
//    return true;
    return newHalfEdges;
}

bool Operator::collapseEdge(HalfEdgeMesh& heMesh, HalfEdge* edge, Eigen::Vector3f newPosition) {
    if (edge == nullptr || edge->twin == nullptr || !checkCollapse(edge)) {
        return false;
    }

    // Update the position of the vertex to the new position
    edge->vertex->position = newPosition;

    // Find all half-edges that need to be updated to point to the merged vertex
    std::vector<HalfEdge*> edgesToUpdate;
    HalfEdge* he = edge->twin;
    do {
        if (he != edge->twin) {
            edgesToUpdate.push_back(he);
        }
        he = he->twin->next;
    } while (he != edge->twin);

    // Now update the half-edges to point from the merged vertex
    for (HalfEdge* heToUpdate : edgesToUpdate) {
        heToUpdate->vertex = edge->vertex;
    }

    // Reassign twin, next, and vertex
    edge->next->twin->twin = edge->next->next->twin;
    edge->next->next->twin->twin = edge->next->twin;
    edge->twin->next->next->twin->twin = edge->twin->next->twin;
    edge->twin->next->twin->twin = edge->twin->next->next->twin;

    edge->vertex->halfEdge = edge->next->next->twin;
    edge->next->next->vertex->halfEdge = edge->next->twin;
    edge->twin->next->next->vertex->halfEdge = edge->twin->next->twin;

    // Remove the half-edges associated with the faces that will be deleted due to the collapse
    HalfEdge* edgeNext = edge->next;
    HalfEdge* edgeTwinNext = edge->twin->next;

    heMesh.halfEdges.remove_if([edgeNext, edgeTwinNext](const HalfEdge& he) {
        return &he == edgeNext || &he == edgeNext->next ||
               &he == edgeTwinNext || &he == edgeTwinNext->next;
    });

    // Remove the vertices associated with the collapse
    Vertex* vertexToRemove = edge->twin->vertex;
    heMesh.vertices.remove_if([vertexToRemove](const Vertex& v) {
        return &v == vertexToRemove;
    });

    // Finally, remove the edge and its twin
    heMesh.halfEdges.remove_if([edge](const HalfEdge& he) {
        return &he == edge || &he == edge->twin;
    });

    return true;
}

bool Operator::checkCollapse(HalfEdge* edge) {
    // Retrieve the vertices at the ends of the edge
    Vertex* v1 = edge->vertex;
    Vertex* v2 = edge->twin->vertex;

    // Use sets to store unique neighboring vertices
    std::set<Vertex*> neighborsOfV1;
    std::set<Vertex*> neighborsOfV2;

    // Helper lambda to add unique neighbors for a given vertex
    auto addNeighbors = [](Vertex* vertex, std::set<Vertex*>& neighbors) {
        HalfEdge* startEdge = vertex->halfEdge;
        HalfEdge* currentEdge = startEdge;
        do {
            if (currentEdge->twin->vertex != vertex) {
                neighbors.insert(currentEdge->twin->vertex);
            }
            currentEdge = currentEdge->twin->next;
        } while (currentEdge != startEdge);
    };

    // Add unique neighbors for both vertices
    addNeighbors(v1, neighborsOfV1);
    addNeighbors(v2, neighborsOfV2);

    // Check for more than two shared neighbors
    std::vector<Vertex*> sharedNeighbors;
    std::set_intersection(neighborsOfV1.begin(), neighborsOfV1.end(),
                          neighborsOfV2.begin(), neighborsOfV2.end(),
                          std::back_inserter(sharedNeighbors));

    if (sharedNeighbors.size() > 2) {
        return false; // Condition 1: More than two shared neighbors
    }

    // Check for a shared neighbor with degree 3
    for (Vertex* sharedNeighbor : sharedNeighbors) {
        int degree = 0;
        HalfEdge* edgeIter = sharedNeighbor->halfEdge;
        do {
            degree++;
            edgeIter = edgeIter->twin->next;
        } while (edgeIter != sharedNeighbor->halfEdge);

        if (degree == 3) {
            return false; // Condition 2: Shared neighbor with degree 3
        }
    }

    return true; // The edge can be collapsed
}

int Operator::vertexDegree(Vertex* vertex) {
    HalfEdge* startEdge = vertex->halfEdge;
    HalfEdge* edge = startEdge;
    int degree = 0;
    do {
        degree++;
        edge = edge->twin->next;
    } while (edge && edge != startEdge);

    return degree;
}

/********** Loop subdivision **********/
void Operator::loopSubdivision(HalfEdgeMesh& heMesh) {
    Operator::SplitResult result = Operator::splitAllEdges(heMesh);
    Operator::flipNewEdges(heMesh, result);
    Operator::updateNewVertices(heMesh, result);
    Operator::updateOldVertices(heMesh, result);
}

Operator::SplitResult Operator::splitAllEdges(HalfEdgeMesh& heMesh) {
    std::vector<HalfEdge*> edgesToSplit;
    std::vector<Vertex*> newVertices;
    std::unordered_set<HalfEdge*> newHalfEdges;

    // Collect all edges to split, including their twins
    for (auto& edge : heMesh.halfEdges) {
        edgesToSplit.push_back(&edge);
        if (edge.twin) {
            edgesToSplit.push_back(edge.twin);
        }
    }

    // Remove duplicate edges from the list, as each edge and its twin would be included twice
    std::sort(edgesToSplit.begin(), edgesToSplit.end());
    edgesToSplit.erase(std::unique(edgesToSplit.begin(), edgesToSplit.end()), edgesToSplit.end());

    // Use a set to keep track of the original edges to prevent splitting newly created edges
    std::unordered_set<HalfEdge*> originalEdges(edgesToSplit.begin(), edgesToSplit.end());

    // Split edges
    for (auto* edge : edgesToSplit) {
        // Ensure only split each original edge and its twin once
        if (originalEdges.find(edge) != originalEdges.end() && originalEdges.find(edge->twin) != originalEdges.end()) {
            size_t prevSize = heMesh.vertices.size();
            std::vector<HalfEdge*> newSplitEdges = splitEdge(heMesh, edge);
            if (heMesh.vertices.size() > prevSize) {
                newVertices.push_back(&heMesh.vertices.back());

                // Remove the edge and its twin from the set after splitting to avoid splitting them again
                originalEdges.erase(edge);
                originalEdges.erase(edge->twin);

                newHalfEdges.insert(newSplitEdges.at(0));
                newHalfEdges.insert(newSplitEdges.at(2));
            }

        }
    }

    Operator::SplitResult result;
    result.newVertices = newVertices;
    result.edgesToFlip = newHalfEdges;

    return result;
}

void Operator::flipNewEdges(HalfEdgeMesh& heMesh, const SplitResult& splitResult) {
    std::vector<Vertex*> newVertices = splitResult.newVertices;
    std::unordered_set<Vertex*> newVerticesSet(newVertices.begin(), newVertices.end());
    for (HalfEdge* edge : splitResult.edgesToFlip) {
        if (newVerticesSet.find(edge->vertex) == newVerticesSet.end() &&
            newVerticesSet.find(edge->twin->vertex) != newVerticesSet.end()){
            flipEdge(heMesh, edge);
        }
    }
}

void Operator::updateNewVertices(HalfEdgeMesh& heMesh, const SplitResult& splitResult) {
    std::vector<Vertex*> newVertices = splitResult.newVertices;
    std::unordered_set<Vertex*> newVerticesSet(newVertices.begin(), newVertices.end());

    std::unordered_set<Vertex*> adjacentVertices;
    std::unordered_set<Vertex*> oppositeVertices;
    for (Vertex* newVertex : newVertices) {
        HalfEdge* startEdge = newVertex->halfEdge;
        HalfEdge* edge = startEdge;
        do {
            if (newVerticesSet.find(edge->twin->vertex) == newVerticesSet.end()) {
                adjacentVertices.insert(edge->twin->vertex);
            }
            else {
                if (newVerticesSet.find(edge->next->twin->next->next->vertex) == newVerticesSet.end() &&
                    newVerticesSet.find(edge->next->twin->next->next->next->vertex) != newVerticesSet.end()) {
                    oppositeVertices.insert(edge->next->twin->next->next->vertex);
                }
            }
            edge = edge->twin->next;
        } while (edge && edge != startEdge);

        newVertex->position = Eigen::Vector3f(0.0, 0.0, 0.0);
        for (Vertex* adjacentVertex : adjacentVertices) {
            newVertex->position += 0.375 * adjacentVertex->position;
        }
        for (Vertex* oppositeVertex : oppositeVertices) {
            newVertex->position += 0.125 * oppositeVertex->position;
        }
        adjacentVertices.clear();
        oppositeVertices.clear();
    }
}

void Operator::updateOldVertices(HalfEdgeMesh& heMesh, const SplitResult& splitResult) {
    std::vector<Vertex*> newVertices = splitResult.newVertices;
    std::unordered_set<Vertex*> newVerticesSet(newVertices.begin(), newVertices.end());
    for (Vertex& v : heMesh.vertices) {
        if (newVerticesSet.find(&v) == newVerticesSet.end()) {
            int neighborCount = 0;
            std::vector<Vertex*> neighborOldVertices;

            HalfEdge* startEdge = v.halfEdge;
            HalfEdge* edge = startEdge;
            do {
                neighborCount++;
                neighborOldVertices.push_back(edge->next->twin->next->twin->next->next->vertex);
                edge = edge->twin->next;
            } while (edge && edge != startEdge);

            float u;
            if (neighborCount == 3) {
                u = 3/16;
            }
            else {
                u = (1.0/neighborCount) * (5.0/8.0 - pow((3.0/8.0 + 1.0/4.0*cos(2*M_PI/neighborCount)), 2));
            }

            v.position = (1 - neighborCount*u) * v.position;
            for (Vertex* neighborOldVertex : neighborOldVertices) {
                v.position += u*neighborOldVertex->position;
            }
        }
    }
}

/********** Quadric Error Simplification **********/
void Operator::quadricErrSimplification(HalfEdgeMesh& heMesh, int targetNum) {
    std::multiset<EdgeErr> edgeQueue;
    for (int i=0; i<targetNum; i++) {
        for (auto& edge : heMesh.halfEdges) {
            Vertex* v1 = edge.vertex;
            Vertex* v2 = edge.twin->vertex;
            Eigen::Matrix4f QSum;
            QSum.setZero();
            QSum = computeQSum(v1) + computeQSum(v2);
            Eigen::Matrix4f QErr = QSum;
            QErr(3, 0) = 0.0f;
            QErr(3, 1) = 0.0f;
            QErr(3, 2) = 0.0f;
            QErr(3, 3) = 1.0f;
            Eigen::Vector3f newPosition = computeMinErrorPoint(QErr);
            float error = computeError(newPosition, QSum);

            EdgeErr edgeErr;
            edgeErr.edgePtr = &edge;
            edgeErr.error = error;
            edgeErr.newPosition = newPosition;
            edgeQueue.insert(edgeErr);
        }

        // Collapse the first edge that can be collapsed
        for (EdgeErr edgeToCollapse : edgeQueue) {
            if (checkCollapse(edgeToCollapse.edgePtr)) {
                collapseEdge(heMesh, edgeToCollapse.edgePtr, edgeToCollapse.newPosition);
                break;
            }
        }

        edgeQueue.clear();
    }
}

Eigen::Matrix4f Operator::computeQ(HalfEdge* edge) {
    HalfEdge* edgeNext = edge->next;
    HalfEdge* edgeNextNext = edge->next->next;

    Eigen::Vector3f v1 = edge->vertex->position;
    Eigen::Vector3f v2 = edgeNext->vertex->position;
    Eigen::Vector3f v3 = edgeNextNext->vertex->position;

    Eigen::Vector3f v1v2 = (v2 -v1).normalized();
    Eigen::Vector3f v1v3 = (v3 -v1).normalized();

    Eigen::Vector3f normal = v1v2.cross(v1v3);
    normal = normal.normalized();

    float d = -(normal.dot(v1));

    Eigen::Vector4f v = Eigen::Vector4f(normal.x(), normal.y(), normal.z(), d);
    Eigen::Matrix4f Q = v * v.transpose();

    return Q;
}

Eigen::Matrix4f Operator::computeQSum(Vertex* vertex) {
    HalfEdge* startEdge = vertex->halfEdge;
    HalfEdge* edge = startEdge;
    Eigen::Matrix4f QSum;
    QSum.setZero();
    do {
        QSum += computeQ(edge);
        edge = edge->twin->next;
    } while (edge && edge != startEdge);

    return QSum;
}

Eigen::Vector3f Operator::computeMinErrorPoint(Eigen::Matrix4f Q) {
    Eigen::Matrix4f Q_pseudo_inverse;
    Eigen::Vector4f pH;
    // Compute the SVD (Singular Value Decomposition) of Q
    Eigen::JacobiSVD<Eigen::MatrixXf> svd(Q, Eigen::ComputeFullU | Eigen::ComputeFullV);
    const Eigen::VectorXf& singularValues = svd.singularValues();
    Eigen::MatrixXf singularValuesInv(Q.rows(), Q.cols());
    singularValuesInv.setZero();

    // Compute the inverse of the singular values matrix
    // Note:Use a threshold for small singular values to avoid division by zero
    const float pinvtoler = 1.e-6; // Pseudo-inverse tolerance
    for (int i = 0; i < singularValues.size(); ++i) {
        if (singularValues(i) > pinvtoler)
            singularValuesInv(i, i) = 1.0f / singularValues(i);
        else
            singularValuesInv(i, i) = 0.0f;
    }

    // Compute the pseudo-inverse as V * S_inv * U^T
    Q_pseudo_inverse = svd.matrixV() * singularValuesInv * svd.matrixU().transpose();

    // Multiply by the vector [0, 0, 0, 1]
    pH = Q_pseudo_inverse * Eigen::Vector4f(0, 0, 0, 1);

    // Convert from homogeneous to 3D coordinates
    Eigen::Vector3f p = Eigen::Vector3f(pH.x(), pH.y(), pH.z());

    return p;
}

float Operator::computeError(Eigen::Vector3f v, Eigen::Matrix4f Q) {
    Eigen::Vector4f vH = Eigen::Vector4f(v.x(), v.y(), v.z(), 1);
    float error = vH.transpose() * Q * vH;
    return error;
}


/********** Isotropic Remeshing **********/
void Operator::isotropicRemeshing(HalfEdgeMesh& heMesh, float weight) {
    Operator::splitLongEdges(heMesh);
    Operator::collapseShortEdges(heMesh);
    Operator::flipForDegree(heMesh);

    for (Vertex& vertex : heMesh.vertices) {
        Operator::tangentialSmoothing(&vertex, weight);
    }
}

void Operator::splitLongEdges(HalfEdgeMesh& heMesh) {
    std::multiset<EdgeLength> edgeLengths;

    float totalLength = 0;
    int counter = 0;
    for (auto& he : heMesh.halfEdges) {
        Eigen::Vector3f v1 = he.vertex->position;
        Eigen::Vector3f v2 = he.twin->vertex->position;
        float length = (v1-v2).norm();

        totalLength += length;
        counter ++;

        EdgeLength edgeLength;
        edgeLength.length = length,
            edgeLength.edgePtr = &he;
        edgeLengths.insert(edgeLength);
    }

    float meanLength = totalLength / counter;

    for (EdgeLength edgeLength : edgeLengths) {
        if (edgeLength.edgePtr == nullptr or edgeLength.edgePtr->twin == nullptr) {
            continue;
        }

        if (edgeLength.length > 1.5*meanLength) {
            std::vector<HalfEdge*> newEdges = splitEdge(heMesh, edgeLength.edgePtr);
        }
    }
}

void Operator::collapseShortEdges(HalfEdgeMesh& heMesh) {
    std::multiset<EdgeLength> edgeLengths;

    float totalLength = 0;
    int counter = 0;
    for (auto& he : heMesh.halfEdges) {
        Eigen::Vector3f v1 = he.vertex->position;
        Eigen::Vector3f v2 = he.twin->vertex->position;
        float length = (v1-v2).norm();

        totalLength += length;
        counter ++;

        EdgeLength edgeLength;
        edgeLength.length = length,
            edgeLength.edgePtr = &he;
        edgeLengths.insert(edgeLength);
    }

    float meanLength = totalLength / counter;

    for (EdgeLength edgeLength : edgeLengths) {
        if (edgeLength.edgePtr == nullptr or edgeLength.edgePtr->twin == nullptr) {
            continue;
        }

        if (edgeLength.length < 0.5*meanLength) {
            if (checkCollapse(edgeLength.edgePtr)) {
                Eigen::Vector3f newPosition = (edgeLength.edgePtr->vertex->position + edgeLength.edgePtr->twin->vertex->position) / 2;
                collapseEdge(heMesh, edgeLength.edgePtr, newPosition);
            }
        }
    }
}

void Operator::flipForDegree(HalfEdgeMesh& heMesh) {
    for (Vertex& vertex : heMesh.vertices) {
        int degree = Operator::vertexDegree(&vertex);

        if (degree > 6) {
            if (checkFlip(vertex.halfEdge)) {
                Operator::flipEdge(heMesh, vertex.halfEdge);
                break;
            }
        }
    }
}


void Operator::tangentialSmoothing(Vertex* vertex, float weight) {
    Eigen::Vector3f neighborTotal = Eigen::Vector3f(0,0,0);
    int counter = 0;
    HalfEdge* startEdge = vertex->halfEdge;
    HalfEdge* edge = startEdge;
    do {
        neighborTotal += edge->twin->vertex->position;
        counter ++;
        edge = edge->twin->next;
    } while (edge && edge != startEdge);

    Eigen::Vector3f v = neighborTotal/counter - vertex->position;
    Eigen::Vector3f n = computeVertexNormal(vertex);
    Eigen::Vector3f differVector = v - n.dot(v) * n;

    vertex->position = vertex->position + weight * differVector;
}

/********** Bilateral Mesh Denoising **********/
void Operator::bilateralMeshDenoising(HalfEdgeMesh& heMesh, float param1, float param2, int rho) {
    for (Vertex& vertex : heMesh.vertices) {
//        std::vector<Eigen::Vector3f> qList;
//        HalfEdge* startEdge = vertex.halfEdge;
//        HalfEdge* edge = startEdge;
//        do {
//            qList.push_back(edge->twin->vertex->position);
//            edge = edge->twin->next;
//        } while (edge && edge != startEdge);
//        int k = qList.size();

        std::vector<Eigen::Vector3f> qList;
        std::queue<std::pair<Vertex*, int>> verticesQueue; // Queue to hold vertices and their depth
        std::unordered_set<Vertex*> visited; // To avoid revisiting vertices

        verticesQueue.push(std::make_pair(&vertex, 0));
        visited.insert(&vertex);

        while (!verticesQueue.empty()) {
            auto front = verticesQueue.front();
            Vertex* currentVertex = front.first;
            int currentDepth = front.second;
            verticesQueue.pop();

            if (currentDepth > 0) {
                qList.push_back(currentVertex->position);
            }

            if (currentDepth < rho) {
                HalfEdge* startEdge = currentVertex->halfEdge;
                HalfEdge* edge = startEdge;
                do {
                    Vertex* nextVertex = edge->twin->vertex;
                    if (visited.find(nextVertex) == visited.end()) {
                        verticesQueue.push(std::make_pair(nextVertex, currentDepth + 1));
                        visited.insert(nextVertex);
                    }
                    edge = edge->twin->next;
                } while (edge && edge != startEdge);
            }
        }
        int k = qList.size();

        Eigen::Vector3f vNormal = computeVertexNormal(&vertex);

        float sum = 0.0;
        float normalizer = 0.0;
        for (int i = 0; i < k; i++) {
            float t = (vertex.position - qList.at(i)).norm();
            float h = vNormal.dot(qList.at(i) - vertex.position);
            const float sigma_c = param1;
            const float sigma_s = param2;
            float wc = exp(-pow(t, 2) / (2 * pow(sigma_c, 2)));
            float ws = exp(-pow(h, 2) / (2 * pow(sigma_s, 2)));
            sum += (wc * ws) * h;
            normalizer += wc * ws;
        }
        vertex.position = vertex.position + vNormal * (sum / normalizer);
    }
}

Eigen::Vector3f Operator::computeFaceNormal(HalfEdge* edge) {
    HalfEdge* edgeNext = edge->next;
    HalfEdge* edgeNextNext = edge->next->next;

    Eigen::Vector3f v1 = edge->vertex->position;
    Eigen::Vector3f v2 = edgeNext->vertex->position;
    Eigen::Vector3f v3 = edgeNextNext->vertex->position;

    Eigen::Vector3f v1v2 = (v2 - v1).normalized();
    Eigen::Vector3f v1v3 = (v3 - v1).normalized();

    Eigen::Vector3f normal = v1v2.cross(v1v3);
    normal.normalize();

    return normal;
}

Eigen::Vector3f Operator::computeVertexNormal(Vertex* vertex) {
    HalfEdge* startEdge = vertex->halfEdge;
    HalfEdge* edge = startEdge;
    Eigen::Vector3f normalSum = Eigen::Vector3f(0,0,0);
    do {
        normalSum += computeFaceNormal(edge);
        edge = edge->twin->next;
    } while (edge && edge != startEdge);

    return normalSum.normalized();
}

float Operator::calculateDistanceToPlane(const Eigen::Vector3f& N, const Eigen::Vector3f& Q, const Eigen::Vector3f& P) {
    // Given a normal vector N and a point Q on the plane
    // Calculate D using point Q (x1, y1, z1)
    float D = -(N.dot(Q));

    // Calculate the distance using point P (x0, y0, z0)
    float distance = std::abs(N.dot(P) + D) / N.norm();

    return distance;
}

void Operator::addNoise(HalfEdgeMesh& heMesh, float magnitude) {
    // Initialize the random number generator once using a non-deterministic seed if possible
    std::random_device rd; // Use std::random_device to generate a seed for more randomness
    std::default_random_engine generator(rd());
    std::uniform_real_distribution<double> distribution(-1.0, 1.0);

    for (Vertex& vertex : heMesh.vertices) {
        Eigen::Vector3f vNormal = computeVertexNormal(&vertex);

        // Generate a random number in the range [-1, 1] using the previously created engine and distribution
        double randomNumber = distribution(generator);
        vertex.position = vertex.position + vNormal * randomNumber * magnitude;
    }
}
