#ifndef SLASUPPORTTREE_HPP
#define SLASUPPORTTREE_HPP

#include <vector>
#include <array>
#include <cstdint>
#include <Eigen/Geometry>

namespace Slic3r {

// Needed types from Point.hpp
typedef int32_t coord_t;
typedef Eigen::Matrix<double,   3, 1, Eigen::DontAlign> Vec3d;
typedef Eigen::Matrix<coord_t,  3, 1, Eigen::DontAlign> Vec3crd;
typedef std::vector<Vec3d>                              Pointf3s;
typedef std::vector<Vec3crd>                            Points3;

class TriangleMesh;

namespace sla {

/// Intermediate struct for a 3D mesh
struct IndexedMesh {
    Pointf3s points;
    Points3  indices;
    std::vector<std::pair<Vec3d, size_t>> weightpoints;
};

/// Generate the 3D support rods for a model intended for SLA print.
void create_support_tree(const IndexedMesh& inpoints, TriangleMesh& output);

}
}

#endif // SLASUPPORTTREE_HPP
