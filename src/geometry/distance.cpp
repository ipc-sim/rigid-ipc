#include "distance.hpp"

#include <ccd/interval.hpp>

namespace ccd {
namespace geometry {

    /// Compute the signed distance between two lines
    template <>
    inline Interval line_line_signed_distance(
        const Eigen::VectorX3<Interval>& line0_point0,
        const Eigen::VectorX3<Interval>& line0_point1,
        const Eigen::VectorX3<Interval>& line1_point0,
        const Eigen::VectorX3<Interval>& line1_point1)
    {
        if (line0_point0.size() != 3) {
            throw NotImplementedError(
                "line_line_signed_distance() not implmeneted in 2D!");
        }
        Eigen::Vector3<Interval> line0_vec = line0_point1 - line0_point0;
        Eigen::Vector3<Interval> line1_vec = line1_point1 - line1_point0;
        Eigen::Vector3<Interval> normal = (line0_vec).cross(line1_vec);
        Interval normal_norm = normal.norm();
        if (boost::numeric::zero_in(normal_norm)) { // parallel lines
            // TODO: Define a signed version of the point line distance
            return point_line_distance(
                line0_point0, line1_point0, line1_point1);
        }
        return (line0_point0 - line1_point0).dot(normal) / normal_norm;
    }

} // namespace geometry
} // namespace ccd
