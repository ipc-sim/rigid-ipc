#include "read_rb_scene.hpp"

#include <unordered_set>

#include <Eigen/Geometry>
#include <boost/filesystem.hpp>
#include <igl/edges.h>
#include <igl/read_triangle_mesh.h>
#include <tbb/parallel_sort.h>

#include <io/read_obj.hpp>
#include <io/serialize_json.hpp>
#include <logger.hpp>
#include <utils/not_implemented_error.hpp>

namespace ccd {
namespace io {

    template <typename T>
    inline bool contains(const std::unordered_set<T>& set, const T& val)
    {
        return set.find(val) != set.end();
    }

    int read_rb_scene_from_str(
        const std::string str, std::vector<physics::RigidBody>& rbs)
    {
        using nlohmann::json;
        json scene = json::parse(str.c_str());
        return read_rb_scene(scene, rbs);
    }

    int read_rb_scene(
        const nlohmann::json& scene, std::vector<physics::RigidBody>& rbs)
    {
        using namespace nlohmann;
        using namespace physics;
        int dim = -1, ndof, angular_dim;
        for (auto& jrb : scene["rigid_bodies"]) {
            // NOTE:
            // All units by default are expressed in standard SI units
            // density: units of Kg/m³ (default density of plastic)
            // position: position of the model origin
            // rotation: degrees as xyz euler angles around the model origin
            // scale: scale the vertices around the model origin
            // linear_velocity: velocity of the body's center
            // angular_velocity: world coordinates in degrees
            // force: applied to center of mass
            // torque: world coordinates in degrees
            json args = R"({
                "mesh": "",
                "vertices": [],
                "edges": [],
                "faces": [],
                "density": 1000.0,
                "is_dof_fixed": [false, false, false, false, false, false],
                "oriented": false,
                "group_id": -1,
                "position": [0.0, 0.0, 0.0],
                "rotation": [0.0, 0.0, 0.0],
                "scale": [1.0, 1.0, 1.0],
                "linear_velocity": [0.0, 0.0, 0.0],
                "angular_velocity": [0.0, 0.0, 0.0],
                "force": [0.0, 0.0, 0.0],
                "torque": [0.0, 0.0, 0.0],
                "enabled": true,
                "type": "dynamic"
            })"_json;
            args.merge_patch(jrb);

            if (!args["enabled"].get<bool>()) {
                continue;
            }

            Eigen::MatrixXd vertices;
            Eigen::MatrixXi faces, edges;

            std::string mesh_fname = args["mesh"].get<std::string>();
            if (mesh_fname != "") {
                boost::filesystem::path mesh_path(mesh_fname);
                if (!exists(mesh_path)) {
                    // TODO: First check a path relative to the input file
                    mesh_path = boost::filesystem::path(__FILE__)
                                    .parent_path() // io
                                    .parent_path() // src
                                    .parent_path() // root
                        / "meshes" / mesh_path;
                }
                spdlog::info("loading mesh: {:s}", mesh_path.string());
                if (mesh_path.extension() == ".obj") {
                    read_obj(mesh_path.string(), vertices, edges, faces);
                } else {
                    igl::read_triangle_mesh(
                        mesh_path.string(), vertices, faces);
                    // Initialize edges
                    if (faces.size()) {
                        igl::edges(faces, edges);
                    }
                }
                assert(faces.size() == 0 || faces.cols() == 3);
            } else {
                // Assumes that edges contains the edges of the faces too.
                from_json<double>(args["vertices"], vertices);
                from_json<int>(args["edges"], edges);
                from_json<int>(args["faces"], faces);
            }

            if (dim == -1) {
                if (vertices.cols() != 0) { // Why would we have an empty body?
                    dim = vertices.cols();
                    ndof = Pose<double>::dim_to_ndof(dim);
                    angular_dim = dim == 2 ? 1 : 3;
                }
            } else if (dim != vertices.cols()) {
                spdlog::error(
                    "Mixing 2D and 3D bodies are not currently allowed.");
                throw NotImplementedError(
                    "Mixing 2D and 3D bodies are not currently allowed.");
            }

            if (dim == 2 && faces.size() != 0) {
                spdlog::warn("Ignoring faces for 2D rigid body.");
                // Delete the faces since they will not be uses
                faces.resize(0, 0);
            }

            Eigen::VectorX3d position;
            from_json<double>(args["position"], position);
            assert(position.size() >= dim);
            position.conservativeResize(dim);

            Eigen::VectorX3d scale;
            // Dimensions overrides a scale argument
            if (args.contains("dimensions")) {
                Eigen::VectorX3d initial_dimensions =
                    (vertices.colwise().maxCoeff()
                     - vertices.colwise().minCoeff())
                        .cwiseAbs();
                for (int i = 0; i < initial_dimensions.size(); i++) {
                    if (initial_dimensions(i) == 0) {
                        initial_dimensions(i) = 1.0;
                    }
                }
                from_json<double>(args["dimensions"], scale);
                assert(scale.size() >= dim);
                scale.conservativeResize(dim);
                scale.array() /= initial_dimensions.array();
            } else if (args["scale"].is_number()) {
                scale.setConstant(dim, args["scale"].get<double>());
            } else {
                from_json<double>(args["scale"], scale);
                assert(scale.size() >= dim);
                scale.conservativeResize(dim);
            }
            vertices = vertices * scale.asDiagonal();

            // Rotate around the models origin NOT the rigid bodies center of
            // mass
            Eigen::VectorX3d rotation;
            from_json<double>(args["rotation"], rotation);
            assert(rotation.size() >= angular_dim);
            rotation.conservativeResize(angular_dim);
            // Convert to radians for easy use later
            rotation *= M_PI / 180.0;
            Eigen::MatrixXX3d R;
            if (rotation.size() == 3) {
                R = (Eigen::AngleAxisd(rotation.z(), Eigen::Vector3d::UnitZ())
                     * Eigen::AngleAxisd(rotation.y(), Eigen::Vector3d::UnitY())
                     * Eigen::AngleAxisd(
                           rotation.x(), Eigen::Vector3d::UnitX()))
                        .toRotationMatrix();
            } else {
                R = Eigen::Rotation2Dd(rotation(0)).toRotationMatrix();
            }
            vertices = vertices * R.transpose();
            rotation.setZero(angular_dim); // Zero initial body rotation

            Eigen::VectorX3d linear_velocity;
            from_json<double>(args["linear_velocity"], linear_velocity);
            assert(linear_velocity.size() >= dim);
            linear_velocity.conservativeResize(dim);

            Eigen::VectorX3d angular_velocity;
            from_json<double>(args["angular_velocity"], angular_velocity);
            assert(angular_velocity.size() >= angular_dim);
            angular_velocity.conservativeResize(angular_dim);
            // Convert to radians for easy use later
            angular_velocity *= M_PI / 180.0;

            Eigen::VectorX3d force;
            from_json<double>(args["force"], force);
            assert(force.size() >= dim);
            force.conservativeResize(dim);

            Eigen::VectorX3d torque;
            from_json<double>(args["torque"], torque);
            assert(torque.size() >= angular_dim);
            torque.conservativeResize(angular_dim);
            // Convert to radians for easy use later
            torque *= M_PI / 180.0;

            Eigen::VectorXb is_dof_fixed;
            if (args["is_dof_fixed"].is_boolean()) {
                is_dof_fixed.setConstant(
                    ndof, args["is_dof_fixed"].get<bool>());
            } else {
                from_json<bool>(args["is_dof_fixed"], is_dof_fixed);
                assert(is_dof_fixed.size() >= ndof);
                is_dof_fixed.conservativeResize(ndof);
            }

            double density = args["density"];
            bool is_oriented = args["oriented"];

            int group_id = args["group_id"];

            RigidBodyType rb_type = args["type"];

            auto rb = physics::RigidBody::from_points(
                vertices, edges, faces, Pose<double>(position, rotation),
                Pose<double>(linear_velocity, angular_velocity),
                Pose<double>(force, torque), density, is_dof_fixed, is_oriented,
                group_id, rb_type);

            rbs.push_back(rb);
        }

        // Adjust the group ids, so the default ones are unique.
        std::unordered_set<int> static_group_ids;
        for (RigidBody& rb : rbs) {
            if (rb.type != RigidBodyType::DYNAMIC) {
                if (rb.group_id >= 0) {
                    static_group_ids.insert(rb.group_id);
                }
                // All static bodies will be given the same group id later
                rb.group_id = -1;
            }
        }
        std::vector<int> taken_ids;
        for (RigidBody& rb : rbs) {
            if (rb.group_id >= 0) {
                taken_ids.push_back(rb.group_id);
            }
        }
        tbb::parallel_sort(taken_ids.begin(), taken_ids.end());
        int taken_id_i = 0;
        int id = 0;
        int static_group_id = -1;
        for (RigidBody& rb : rbs) {
            bool in_static_group = rb.type != RigidBodyType::DYNAMIC
                || contains(static_group_ids, rb.group_id);
            if (static_group_id >= 0 && in_static_group) {
                // All static bodies are in the same group
                rb.group_id = static_group_id;
            } else if (rb.group_id < 0) {
                // Find the next free id
                while (taken_id_i < taken_ids.size()
                       && id >= taken_ids[taken_id_i]) {
                    if (id == taken_ids[taken_id_i]) {
                        id++;
                    }
                    taken_id_i++;
                }
                rb.group_id = id++;
                assert(!std::binary_search(
                    taken_ids.begin(), taken_ids.end(), rb.group_id));
            }

            assert(rb.group_id >= 0);
            if (static_group_id < 0 && in_static_group) {
                // All static bodies are in the same group
                static_group_id = rb.group_id;
            }
        }

        return dim;
    }

} // namespace io
} // namespace ccd
