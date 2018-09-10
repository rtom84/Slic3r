#ifndef SLASUPPORTTREE_HPP
#define SLASUPPORTTREE_HPP

#include <vector>
#include <Eigen/Geometry>

namespace Slic3r {

// Needed types from Point.hpp
typedef Eigen::Matrix<double,   3, 1, Eigen::DontAlign> Vec3d;
typedef std::vector<Vec3d>                              Pointf3s;

class TriangleMesh;

namespace sla {

class Vec3dN: public Vec3d {
    Vec3d n_;
public:

    template<class...Args> Vec3dN(const Vec3d& normal, Args&&...args):
        Vec3d(std::forward<Args>(args)...), n_(normal) {}

    const Vec3d& normal() const { return n_; }
};

using WeightPoints = std::vector<Vec3dN>;

/// Generate the 3D support rods for a model intended for SLA print.
void create_support_tree(const WeightPoints& inpoints, TriangleMesh& output);

}
}

#endif // SLASUPPORTTREE_HPP
