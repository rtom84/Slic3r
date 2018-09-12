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

    // A vector of pairs: first is the point of interest regarding the support
    // generation, second is an index into "indices" that corresponds to the
    // triangle which owns the POI.
    // The POI can be in baricentric or world coordinates... we should decide
    std::vector<std::pair<Vec3d, size_t>> weightpoints;
};

struct SupportConfig {
    // Radius in radians of the pointing side of the head.
    double head_front_radius_rad;

    // Radius of the back side of the 3d arrow.
    double head_back_radius_rad;

    // Width in mm from the back sphere center to the front sphere center.
    double head_width_mm;

    // Radius in radians of the support pillars.
    double pillar_radius_rad;

    // Radius in radians of the pillar base.
    double base_radius_rad;

    // The height of the pillar base cone in mm.
    double base_height_mm;
};

/// Generate the 3D support rods for a model intended for SLA print.
void create_support_tree(const IndexedMesh& inpoints,
                         TriangleMesh& output,
                         const SupportConfig& cfg);

}
}

#endif // SLASUPPORTTREE_HPP
