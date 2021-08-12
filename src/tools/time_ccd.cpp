#include <CLI/CLI.hpp>

#include <ghc/fs_std.hpp> // filesystem
#include <igl/Timer.h>
#include <igl/edges.h>

#include <ccd/ccd.hpp>
#include <ccd/piecewise_linear/time_of_impact.hpp>
#include <ccd/redon/time_of_impact.hpp>
#include <ccd/rigid/time_of_impact.hpp>
#include <io/serialize_json.hpp>

using namespace ipc;
using namespace ipc::rigid;

RigidBody create_body(
    const Eigen::MatrixXd& vertices,
    const Eigen::MatrixXi& edges,
    const Eigen::MatrixXi& faces)
{
    static int id = 0;
    int dim = vertices.cols();
    PoseD pose = PoseD::Zero(dim);
    RigidBody rb(
        vertices, edges, faces, pose,
        /*velocity=*/PoseD::Zero(pose.dim()),
        /*force=*/PoseD::Zero(pose.dim()),
        /*denisty=*/1.0,
        /*is_dof_fixed=*/VectorMax6b::Zero(pose.ndof()),
        /*oriented=*/false,
        /*group_id=*/id++);
    rb.vertices = vertices; // Cancel out the inertial rotation for testing
    rb.pose.position.setZero();
    rb.pose.rotation.setZero();
    return rb;
}

RigidBody
create_body(const Eigen::MatrixXd& vertices, const Eigen::MatrixXi& edges)
{
    return create_body(vertices, edges, Eigen::MatrixXi());
}

int main(int argc, char* argv[])
{
    set_logger_level(static_cast<spdlog::level::level_enum>(6));
    igl::Timer timer;

    std::ofstream fv_csv("fv.csv");
    fv_csv << "redon,rigid,pl\n";
    std::ofstream ee_csv("ee.csv");
    ee_csv << "redon,rigid,pl\n";

    std::array<TrajectoryType, 3> traj_types = { { REDON, RIGID,
                                                   PIECEWISE_LINEAR } };
    static double REDON_TOL = 1e-4, RIGID_TOL = 1e-4;

    for (int i = 0; i < 60; i++) {
        fs::path data_path = fs::path(__FILE__).parent_path() / "ccd-queries"
            / fmt::format("ccd-queries-{:03d}.json", i);

        std::ifstream input(data_path.string());
        nlohmann::json queries_json = nlohmann::json::parse(input);

        std::vector<nlohmann::json> queries = queries_json["queries"];

        int counter = 0;
        for (const auto& query : queries) {
            std::cout << counter++ << "\r" << std::flush;
            // if (counter > 10) {
            //     break;
            // }

            Eigen::MatrixXd bodyA_vertices, bodyB_vertices;
            Eigen::MatrixXi bodyA_edges, bodyB_edges, bodyA_faces, bodyB_faces;
            PoseD bodyA_pose_t0, bodyA_pose_t1;
            PoseD bodyB_pose_t0, bodyB_pose_t1;

            std::string ccd_type = query["type"];

            nlohmann::json poseA_t0_json, poseA_t1_json;
            nlohmann::json poseB_t0_json, poseB_t1_json;

            if (ccd_type == "ee") {
                Eigen::VectorXd tmp;

                bodyA_vertices.resize(2, 3);
                from_json(query["edge0"]["vertex0"], tmp);
                bodyA_vertices.row(0) = tmp;
                from_json(query["edge0"]["vertex1"], tmp);
                bodyA_vertices.row(1) = tmp;

                bodyB_vertices.resize(2, 3);
                from_json(query["edge1"]["vertex0"], tmp);
                bodyB_vertices.row(0) = tmp;
                from_json(query["edge1"]["vertex1"], tmp);
                bodyB_vertices.row(1) = tmp;

                bodyA_edges.resize(1, 2);
                bodyA_edges.row(0) << 0, 1;
                bodyB_edges.resize(1, 2);
                bodyB_edges.row(0) << 0, 1;

                poseA_t0_json = query["edge0"]["pose_t0"];
                poseA_t1_json = query["edge0"]["pose_t1"];
                poseB_t0_json = query["edge1"]["pose_t0"];
                poseB_t1_json = query["edge1"]["pose_t1"];
            } else if (ccd_type == "fv") {
                Eigen::VectorXd tmp;

                bodyA_vertices.resize(3, 3);
                from_json(query["face"]["vertex0"], tmp);
                bodyA_vertices.row(0) = tmp;
                from_json(query["face"]["vertex1"], tmp);
                bodyA_vertices.row(1) = tmp;
                from_json(query["face"]["vertex2"], tmp);
                bodyA_vertices.row(2) = tmp;

                bodyB_vertices.resize(1, 3);
                from_json(query["vertex"]["vertex"], tmp);
                bodyB_vertices.row(0) = tmp;

                bodyA_faces.resize(1, 3);
                bodyA_faces.row(0) << 0, 1, 2;
                bodyA_edges.resize(3, 2);
                bodyA_edges.row(0) << 0, 1;
                bodyA_edges.row(1) << 1, 2;
                bodyA_edges.row(2) << 2, 0;

                poseA_t0_json = query["face"]["pose_t0"];
                poseA_t1_json = query["face"]["pose_t1"];
                poseB_t0_json = query["vertex"]["pose_t0"];
                poseB_t1_json = query["vertex"]["pose_t1"];
            }

            from_json(poseA_t0_json["position"], bodyA_pose_t0.position);
            from_json(poseA_t0_json["rotation"], bodyA_pose_t0.rotation);
            from_json(poseA_t1_json["position"], bodyA_pose_t1.position);
            from_json(poseA_t1_json["rotation"], bodyA_pose_t1.rotation);
            from_json(poseB_t0_json["position"], bodyB_pose_t0.position);
            from_json(poseB_t0_json["rotation"], bodyB_pose_t0.rotation);
            from_json(poseB_t1_json["position"], bodyB_pose_t1.position);
            from_json(poseB_t1_json["rotation"], bodyB_pose_t1.rotation);

            RigidBody bodyA =
                create_body(bodyA_vertices, bodyA_edges, bodyA_faces);
            RigidBody bodyB =
                create_body(bodyB_vertices, bodyB_edges, bodyB_faces);

            double toi = 1;
            bool is_impacting = false;

            std::array<double, 3> timings;
            for (int i = 0; i < traj_types.size(); i++) {

                timer.start();
                if (ccd_type == "ee") {
                    switch (traj_types[i]) {
                    case REDON:
                        compute_edge_edge_time_of_impact_redon(
                            bodyA, bodyA_pose_t0, bodyA_pose_t1,
                            /*edgeA_id=*/0, //
                            bodyB, bodyB_pose_t0, bodyB_pose_t1,
                            /*edgeB_id=*/0, //
                            toi, /*double earliest_toi=*/1, REDON_TOL);
                    case RIGID:
                        compute_edge_edge_time_of_impact(
                            bodyA, bodyA_pose_t0, bodyA_pose_t1,
                            /*edgeA_id=*/0, //
                            bodyB, bodyB_pose_t0, bodyB_pose_t1,
                            /*edgeB_id=*/0, //
                            toi, /*double earliest_toi=*/1, RIGID_TOL);
                    case PIECEWISE_LINEAR:
                        compute_piecewise_linear_edge_edge_time_of_impact(
                            bodyA, bodyA_pose_t0, bodyA_pose_t1,
                            /*edgeA_id=*/0, //
                            bodyB, bodyB_pose_t0, bodyB_pose_t1,
                            /*edgeB_id=*/0, //
                            toi);
                        break;
                    default:
                        break;
                    }
                } else if (ccd_type == "fv") {
                    switch (traj_types[i]) {
                    case REDON:
                        compute_face_vertex_time_of_impact_redon(
                            bodyB, bodyB_pose_t0, bodyB_pose_t1,
                            /*vertex_id=*/0, // Vertex body
                            bodyA, bodyA_pose_t0, bodyA_pose_t1,
                            /*face_id=*/0, // Face body
                            toi, /*double earliest_toi=*/1, REDON_TOL);
                    case RIGID:
                        compute_face_vertex_time_of_impact(
                            bodyB, bodyB_pose_t0, bodyB_pose_t1,
                            /*vertex_id=*/0, // Vertex body
                            bodyA, bodyA_pose_t0, bodyA_pose_t1,
                            /*face_id=*/0, // Face body
                            toi, /*double earliest_toi=*/1, RIGID_TOL);
                    case PIECEWISE_LINEAR:
                        compute_piecewise_linear_face_vertex_time_of_impact(
                            bodyB, bodyB_pose_t0, bodyB_pose_t1,
                            /*vertex_id=*/0, // Vertex body
                            bodyA, bodyA_pose_t0, bodyA_pose_t1,
                            /*face_id=*/0, // Face body
                            toi);
                        break;
                    default:
                        break;
                    }
                }
                timer.stop();

                timings[i] = timer.getElapsedTimeInMicroSec();
            }

            std::string line = fmt::format(
                "{:g},{:g},{:g}\n", timings[0], timings[1], timings[2]);
            if (ccd_type == "ee") {
                ee_csv << line;
            } else if (ccd_type == "fv") {
                fv_csv << line;
            }
        }

        fv_csv.flush();
        ee_csv.flush();
    }
}
