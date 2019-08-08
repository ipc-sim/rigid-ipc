#include "read_rb_scene.hpp"

#include <fstream>
#include <io/serialize_json.hpp>
#include <iostream>

namespace ccd {
namespace io {

    void read_rb_scene_from_str(
        const std::string str, std::vector<physics::RigidBody>& rbs)
    {
        using nlohmann::json;
        json scene = json::parse(str.c_str());
        return read_rb_scene(scene, rbs);
    }

    void read_rb_scene(
        const nlohmann::json& scene, std::vector<physics::RigidBody>& rbs)
    {
        using namespace nlohmann;
        for (auto& jrb : scene["rigid_bodies"]) {
            json args = R"({
                  "vertices":[],
                  "edges":[],
                  "velocity":[0.0,0.0,0.0],
                  "position":[0.0,0.0],
                  "theta":0.0,
                  "is_dof_fixed":[false,false,false]
                  })"_json;
            args.merge_patch(jrb);

            Eigen::MatrixXd vertices;
            from_json(args["vertices"], vertices);

            Eigen::VectorXd xy_position;
            from_json(args["position"], xy_position);

            Eigen::MatrixXi edges;
            from_json(args["edges"], edges);

            Eigen::VectorXd velocity;
            from_json(args["velocity"], velocity);

            Eigen::VectorXb is_dof_fixed;
            from_json(args["is_dof_fixed"], is_dof_fixed);

            double theta = args["theta"].get<double>() * M_PI / 180.0;
            Eigen::Vector3d position;
            position.head(2) = xy_position;
            position(2) = theta;
            auto rb = physics::RigidBody::from_velocity(
                vertices, edges, position, velocity, is_dof_fixed);



            rbs.push_back(rb);
        }
    }

} // namespace io
} // namespace ccd
