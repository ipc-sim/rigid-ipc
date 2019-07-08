#include "read_rb_scene.hpp"

#include <fstream>
#include <iostream>
#include <io/serialize_json.hpp>

namespace ccd {
namespace io {

    bool is_rb_scene(const nlohmann::json& scene)
    {
        auto scene_type = scene["scene_type"];
        return scene_type.is_string()
            && scene_type.get<std::string>().compare("rigid-body") == 0;
    }

    void read_rb_scene(
        const std::string filename, std::vector<physics::RigidBody>& rbs)
    {
        using nlohmann::json;

        std::ifstream input(filename);
        json scene = json::parse(input);
        return read_rb_scene(scene, rbs);
    }

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
        assert(std::string("rigid-body").compare(scene["scene_type"]) == 0);
        for (auto& jrb : scene["rigid_bodies"]) {
            Eigen::MatrixXd vertices;
            from_json(jrb["vertices"], vertices);
            Eigen::MatrixXi edges;
            from_json(jrb["edges"], edges);
            Eigen::VectorXd velocity;
            from_json(jrb["velocity"], velocity);

            auto rb = physics::RigidBody::Centered(vertices, edges, velocity);
            rbs.push_back(rb);
        }
    }

} // namespace io
} // namespace ccd
