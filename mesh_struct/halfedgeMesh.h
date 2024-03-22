#ifndef HALFEDGE_H
#define HALFEDGE_H

#include <list>
#include <Eigen/Core>
#include <unordered_map>
#include <unordered_set>

struct HalfEdge;
struct Vertex;

struct HalfEdge {
    HalfEdge* twin;
    HalfEdge* next;
    Vertex* vertex;
    // Adding a constructor for convenience
    HalfEdge() : twin(nullptr), next(nullptr), vertex(nullptr) {}
};

struct Vertex {
    Eigen::Vector3f position;
    HalfEdge* halfEdge;
    // Adding a constructor for convenience
    Vertex(const Eigen::Vector3f& pos) : position(pos), halfEdge(nullptr) {}
};

struct pair_hash {
    template <class T1, class T2>
    std::size_t operator () (const std::pair<T1, T2>& pair) const {
        return std::hash<T1>()(pair.first) ^ std::hash<T2>()(pair.second);
    }
};

class HalfEdgeMesh {
public:
    HalfEdgeMesh(const std::vector<Eigen::Vector3f>& vertices,
                 const std::vector<Eigen::Vector3i>& faces);
    ~HalfEdgeMesh();

    std::list<Vertex> vertices;
    std::list<HalfEdge> halfEdges;

    void convertToMeshFormat(std::vector<Eigen::Vector3f>& outVertices,
                             std::vector<Eigen::Vector3i>& outFaces);

    bool validate() const;

private:
    void buildHalfEdgeStructure(const std::vector<Eigen::Vector3f>& _vertices,
                                const std::vector<Eigen::Vector3i>& _faces);
};

#endif // HALFEDGE_H
