#include "read_rb_scene.hpp"

#include <unordered_set>

#include <Eigen/Geometry>
#include <ghc/fs_std.hpp> // filesystem
#include <igl/edges.h>
#include <igl/facet_components.h>
#include <igl/PI.h>
#include <igl/read_triangle_mesh.h>
#include <igl/remove_unreferenced.h>
#include <tbb/parallel_sort.h>

#include <io/read_obj.hpp>
#include <io/serialize_json.hpp>
#include <logger.hpp>
#include <utils/not_implemented_error.hpp>

namespace ipc::rigid {

template <typename T>
inline bool contains(const std::unordered_set<T>& set, const T& val)
{
    return set.find(val) != set.end();
}

bool read_rb_scene_from_str(const std::string str, std::vector<RigidBody>& rbs)
{
    using nlohmann::json;
    json scene = json::parse(str.c_str());
    return read_rb_scene(scene, rbs);
}

VectorMax3d read_angular_field(const nlohmann::json& field, int dim)
{
    int angular_dim = dim == 2 ? 1 : 3;

    VectorMax3d v;
    if (field.is_number()) {
        assert(dim == 2);
        v.setZero(angular_dim);
        v[0] = field.get<double>();
    } else {
        assert(field.is_array());
        from_json(field, v);
        assert(v.size() >= angular_dim);
        v.conservativeResize(angular_dim);
    }

    // Convert to radians for easy use later
    v *= igl::PI / 180.0;

    return v;
}

bool read_rb_scene(const nlohmann::json& scene, std::vector<RigidBody>& rbs)
{
    using namespace nlohmann;
    int dim = -1, ndof, angular_dim;

    std::unordered_map<std::string, int> rb_name_to_count;

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
                "type": "dynamic",
                "kinematic_max_time": -1,
                "kinematic_poses": [],
                "split_components": false
            })"_json;
        args.merge_patch(jrb);
        if (args["kinematic_max_time"] < 0) {
            args["kinematic_max_time"] =
                std::numeric_limits<double>::infinity();
        }

        if (!args["enabled"].get<bool>()) {
            continue;
        }

        Eigen::MatrixXd vertices;
        Eigen::MatrixXi faces, edges;
        std::string rb_name;

        std::string mesh_fname = args["mesh"].get<std::string>();
        if (mesh_fname != "") {
            fs::path mesh_path(mesh_fname);
            if (!exists(mesh_path)) {
                // TODO: First check a path relative to the input file
                mesh_path = fs::path(RIGID_IPC_MESHES_DIR) / mesh_path;
            }
            spdlog::info("loading mesh: {:s}", mesh_path.string());
            bool success;
            if (mesh_path.extension() == ".obj") {
                success = read_obj(mesh_path.string(), vertices, edges, faces);
            } else {
                success = igl::read_triangle_mesh(
                    mesh_path.string(), vertices, faces);
                // Initialize edges
                if (faces.size()) {
                    igl::edges(faces, edges);
                }
            }
            assert(faces.size() == 0 || faces.cols() == 3);
            if (!success) {
                return false;
            }

            rb_name = mesh_path.stem().string();
        } else {
            // Assumes that edges contains the edges of the faces too.
            from_json(args["vertices"], vertices);
            from_json(args["edges"], edges);
            from_json(args["faces"], faces);
            rb_name = "RigidBody";
        }

        if (dim == -1) {
            if (vertices.cols() != 0) { // Why would we have an empty body?
                dim = vertices.cols();
                ndof = PoseD::dim_to_ndof(dim);
                angular_dim = dim == 2 ? 1 : 3;
            }
        } else if (dim != vertices.cols()) {
            spdlog::error("Mixing 2D and 3D bodies are not currently allowed.");
            throw NotImplementedError(
                "Mixing 2D and 3D bodies are not currently allowed.");
        }

        if (dim == 2 && faces.size() != 0) {
            spdlog::warn("Ignoring faces for 2D rigid body.");
            // Delete the faces since they will not be uses
            faces.resize(0, 0);
        }

        VectorMax3d position;
        from_json(args["position"], position);
        assert(position.size() >= dim);
        position.conservativeResize(dim);

        VectorMax3d scale;
        // Dimensions overrides a scale argument
        if (args.contains("dimensions")) {
            VectorMax3d initial_dimensions =
                (vertices.colwise().maxCoeff() - vertices.colwise().minCoeff())
                    .cwiseAbs();
            initial_dimensions =
                (initial_dimensions.array() == 0).select(1, initial_dimensions);
            from_json(args["dimensions"], scale);
            assert(scale.size() >= dim);
            scale.conservativeResize(dim);
            scale.array() /= initial_dimensions.array();
        } else if (args["scale"].is_number()) {
            scale.setConstant(dim, args["scale"].get<double>());
        } else {
            from_json(args["scale"], scale);
            assert(scale.size() >= dim);
            scale.conservativeResize(dim);
        }
        vertices *= scale.asDiagonal();

        // Rotate around the models origin NOT the rigid bodies center of mass
        VectorMax3d rotation = read_angular_field(args["rotation"], dim);
        MatrixMax3d R;
        if (rotation.size() == 3) {
            R = (Eigen::AngleAxisd(rotation.z(), Eigen::Vector3d::UnitZ())
                 * Eigen::AngleAxisd(rotation.y(), Eigen::Vector3d::UnitY())
                 * Eigen::AngleAxisd(rotation.x(), Eigen::Vector3d::UnitX()))
                    .toRotationMatrix();
        } else {
            R = Eigen::Rotation2Dd(rotation(0)).toRotationMatrix();
        }
        vertices = vertices * R.transpose();
        rotation.setZero(angular_dim); // Zero initial body rotation

        VectorMax3d linear_velocity;
        from_json(args["linear_velocity"], linear_velocity);
        assert(linear_velocity.size() >= dim);
        linear_velocity.conservativeResize(dim);

        VectorMax3d angular_velocity =
            read_angular_field(args["angular_velocity"], dim);

        VectorMax3d force;
        from_json(args["force"], force);
        assert(force.size() >= dim);
        force.conservativeResize(dim);

        VectorMax3d torque = read_angular_field(args["torque"], dim);

        VectorXb is_dof_fixed;
        if (args["is_dof_fixed"].is_boolean()) {
            is_dof_fixed.setConstant(ndof, args["is_dof_fixed"].get<bool>());
        } else {
            from_json(args["is_dof_fixed"], is_dof_fixed);
            assert(is_dof_fixed.size() >= ndof);
            is_dof_fixed.conservativeResize(ndof);
        }

        double density = args["density"];
        bool is_oriented = args["oriented"];

        int group_id = args["group_id"];

        RigidBodyType rb_type = args["type"];
        double kinematic_max_time = args["kinematic_max_time"];

        std::vector<json> json_kinematic_poses = args["kinematic_poses"];
        std::deque<PoseD> kinematic_poses;
        for (const auto& json_pose : json_kinematic_poses) {
            PoseD pose = PoseD::Zero(dim);
            if (json_pose.contains("position")) {
                from_json(json_pose["position"], pose.position);
            }
            if (json_pose.contains("rotation")) {
                from_json(json_pose["rotation"], pose.rotation);
            }
            kinematic_poses.push_back(pose);
        }

        if (args["split_components"].get<bool>()) {
            // TODO: Handle codimensional edges too
            assert(faces.cols() == 3);
            Eigen::VectorXi C;
            igl::facet_components(faces, C);
            int num_components = C.maxCoeff();
            std::vector<std::vector<int>> CFs(num_components + 1);
            for (int j = 0; j < faces.cols(); j++) {
                for (int i = 0; i < faces.rows(); i++) {
                    CFs[C[i]].push_back(faces(i, j));
                }
            }

            for (int ci = 0; ci < CFs.size(); ci++) {
                Eigen::MatrixXi F = Eigen::Map<Eigen::MatrixXi>(
                    CFs[ci].data(), CFs[ci].size() / 3, 3);
                Eigen::MatrixXd CV;
                Eigen::MatrixXi CF;
                Eigen::VectorXi I;
                igl::remove_unreferenced(vertices, F, CV, CF, I);
                Eigen::MatrixXi CE;
                igl::edges(CF, CE);
                // WARNING: angular velocity and torque will be around the
                // components center of mass not the entire meshes.
                rbs.emplace_back(
                    CV, CE, CF, PoseD(position, VectorMax3d::Zero(angular_dim)),
                    PoseD(linear_velocity, angular_velocity),
                    PoseD(force, torque), density, is_dof_fixed, is_oriented,
                    group_id, rb_type, kinematic_max_time, kinematic_poses);
                rbs.back().name = fmt::format("{}-part{:03d}", rb_name, ci);
            }
        } else {
            rbs.emplace_back(
                vertices, edges, faces,
                PoseD(position, VectorMax3d::Zero(angular_dim)),
                PoseD(linear_velocity, angular_velocity), PoseD(force, torque),
                density, is_dof_fixed, is_oriented, group_id, rb_type,
                kinematic_max_time, kinematic_poses);
            rbs.back().name = rb_name;
        }
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

    return true;
}

} // namespace ipc::rigid
