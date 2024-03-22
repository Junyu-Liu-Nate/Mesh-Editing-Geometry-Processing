#pragma once

#include "mesh_struct/halfedgeMesh.h"

class Operator
{
public:
    Operator();

    //TODO: Implement basic operators:
    //      - edge flip, edge split, edge collapse
    static HalfEdge* findFlipableEdge(HalfEdgeMesh& mesh);
    static bool flipEdge(HalfEdgeMesh& heMesh, HalfEdge* edge);
    static std::pair<HalfEdge*, HalfEdge*> flipEdgeReturn(HalfEdgeMesh& heMesh, HalfEdge* edge);
    static bool checkFlip(HalfEdge* edge);
//    static bool splitEdge(HalfEdgeMesh& heMesh, HalfEdge* edge);
    static std::vector<HalfEdge*> splitEdge(HalfEdgeMesh& heMesh, HalfEdge* edge);
    static bool collapseEdge(HalfEdgeMesh& heMesh, HalfEdge* edge, Eigen::Vector3f newPosition);
    static bool checkCollapse(HalfEdge* edge);

    static int vertexDegree(Vertex* vertex);

    //TODO: Implement more complex opreations
    //      Use the basic operators
    /********** Loop subdivision **********/
    static void loopSubdivision(HalfEdgeMesh& heMesh);
    struct SplitResult {
        std::vector<Vertex*> newVertices;
        std::unordered_set<HalfEdge*> edgesToFlip;
    };
    static SplitResult splitAllEdges(HalfEdgeMesh& heMesh);
    static void flipNewEdges(HalfEdgeMesh& heMesh, const SplitResult& splitResult);
    static void updateNewVertices(HalfEdgeMesh& heMesh, const SplitResult& splitResult);
    static void updateOldVertices(HalfEdgeMesh& heMesh, const SplitResult& splitResult);

    /********** Quadric Error Simplification **********/
    static void quadricErrSimplification(HalfEdgeMesh& heMesh, int targetNum);
    struct VertexQ {
        Eigen::Matrix4f Q;
        Vertex* vertexPtr;
    };
    struct EdgeErr {
        float error;
        Eigen::Vector3f newPosition;
        HalfEdge* edgePtr;

        // The less-than operator is required for ordering in the multiset
        bool operator<(const EdgeErr& other) const {
            return error < other.error;
        }
    };
    static Eigen::Matrix4f computeQ(HalfEdge* edge);
    static Eigen::Matrix4f computeQSum(Vertex* v);
    static Eigen::Vector3f computeMinErrorPoint(Eigen::Matrix4f Q);
    static float computeError(Eigen::Vector3f v, Eigen::Matrix4f Q);

    /********** Isotropic Remeshing **********/
    static void isotropicRemeshing(HalfEdgeMesh& heMesh, float weight);
    struct EdgeLength {
        float length;
        HalfEdge* edgePtr;

        bool operator<(const EdgeLength& other) const {
            return length < other.length;
        }
    };
    static void splitLongEdges(HalfEdgeMesh& heMesh);
    static void collapseShortEdges(HalfEdgeMesh& heMesh);
    static void flipForDegree(HalfEdgeMesh& heMesh);
    static void tangentialSmoothing(Vertex* vertex, float weight);

    /********** Bilateral Mesh Denoising **********/
    static void bilateralMeshDenoising(HalfEdgeMesh& heMesh, float param1, float param2, int rho);
    static Eigen::Vector3f computeFaceNormal(HalfEdge* edge);
    static Eigen::Vector3f computeVertexNormal(Vertex* vertex);
    static float calculateDistanceToPlane(const Eigen::Vector3f& N, const Eigen::Vector3f& Q, const Eigen::Vector3f& P);

    static void addNoise(HalfEdgeMesh& heMesh, float magnitude);
private:
};
