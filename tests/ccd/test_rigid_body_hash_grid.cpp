#include <catch2/catch.hpp>

#include <filesystem>
namespace fs = std::filesystem;

#include <igl/Timer.h>
#include <nlohmann/json.hpp>

#include <io/serialize_json.hpp>
#include <logger.hpp>
#include <physics/pose.hpp>

using namespace ipc;
using namespace ipc::rigid;

TEST_CASE("2D rigid body hash grid", "[hashgrid][rigid_body][2D]")
{
    // TODO
}

TEST_CASE("3D rigid body hash grid", "[hashgrid][rigid_body][3D]")
{
    // TODO
}

void compute_scene_conservative_bbox(
    const std::vector<nlohmann::json>& bodies,
    Eigen::Vector3d& scene_min,
    Eigen::Vector3d& scene_max)
{
    scene_min.setConstant(std::numeric_limits<double>::infinity());
    scene_max.setConstant(-std::numeric_limits<double>::infinity());

    for (const auto& body : bodies) {
        Eigen::MatrixXd V;
        from_json(body["vertices"], V);
        double radius = V.rowwise().norm().maxCoeff();

        Eigen::Vector3d p0, p1;
        Eigen::Vector3d r0, r1;
        from_json(body["pose_t0"]["position"], p0);
        from_json(body["pose_t0"]["rotation"], r0);
        from_json(body["pose_t1"]["position"], p1);
        from_json(body["pose_t1"]["rotation"], r1);

        if ((r0.array() == r1.array()).all()) {
            auto R = construct_rotation_matrix(VectorMax3d(r0));
            Eigen::MatrixXd V0 = (V * R.transpose()).rowwise() + p0.transpose();
            scene_min = scene_min.cwiseMin(V0.colwise().minCoeff().transpose());
            scene_max = scene_max.cwiseMax(V0.colwise().maxCoeff().transpose());
            Eigen::MatrixXd V1 = (V * R.transpose()).rowwise() + p1.transpose();
            scene_min = scene_min.cwiseMin(V1.colwise().minCoeff().transpose());
            scene_max = scene_max.cwiseMax(V1.colwise().maxCoeff().transpose());
        } else {
            scene_min = scene_min.cwiseMin((p0.array() - radius).matrix());
            scene_max = scene_max.cwiseMax((p0.array() + radius).matrix());
            scene_min = scene_min.cwiseMin((p1.array() - radius).matrix());
            scene_max = scene_max.cwiseMax((p1.array() + radius).matrix());
        }
    }
}

Vector3I compute_scene_bbox(
    const std::vector<nlohmann::json>& bodies, int num_subdivisions)
{
    Vector3I scene_bbox = Vector3I::Constant(Interval::empty());

    double ti0 = 0;
    for (int i = 0; i < num_subdivisions; i++) {
        double ti1 = ti0 + 1.0 / num_subdivisions;
        Interval t(ti0, ti1);

        for (const auto& body : bodies) {
            Eigen::Vector3d p0d, p1d;
            Eigen::Vector3d r0d, r1d;
            from_json(body["pose_t0"]["position"], p0d);
            from_json(body["pose_t0"]["rotation"], r0d);
            from_json(body["pose_t1"]["position"], p1d);
            from_json(body["pose_t1"]["rotation"], r1d);

            Vector3I p0 = p0d.cast<Interval>(), p1 = p1d.cast<Interval>();
            Vector3I r0 = r0d.cast<Interval>(), r1 = r1d.cast<Interval>();

            Vector3I p = (p1 - p0) * t + p0;
            Vector3I r = (r1 - r0) * t + r0;
            auto R = construct_rotation_matrix(VectorMax3I(r));

            Eigen::MatrixXd V;
            from_json(body["vertices"], V);
            MatrixXI VI = (V * R.transpose()).rowwise() + p.transpose();

            for (int i = 0; i < VI.rows(); i++) {
                for (int j = 0; j < VI.cols(); j++) {
                    scene_bbox(j) = hull(scene_bbox(j), VI(i, j));
                }
            }
        }

        ti0 = ti1;
    }

    return scene_bbox;
}

// TEST_CASE("Large hashgrid", "[hashgrid][rigid_body][3D]")
// {
//     // Load json file
//     std::string filename = "large-rb-hashgrid-001.json";
//     fs::path data_path =
//         fs::path(__FILE__).parent_path().parent_path() /
//         "data" / "large-rb-hashgrid" / filename;
//     std::ifstream input(data_path.string());
//     REQUIRE(input.is_open());
//     nlohmann::json json = nlohmann::json::parse(input);
//
//     std::vector<nlohmann::json> bodies = json["bodies"];
//
//     Eigen::Vector3d conservative_min, conservative_max;
//     compute_scene_conservative_bbox(bodies, conservative_min,
//     conservative_max); conservative_min.array() -= 1e-8;
//     conservative_max.array() += 1e-8;
//
//     fmt::print("num_subdivisions,scene_bbox_diag,time\n");
//     igl::Timer timer;
//     for (int n = 1; n < 20; n++) {
//         timer.start();
//         Vector3I scene_bbox = compute_scene_bbox(bodies, n);
//         timer.stop();
//         Eigen::Vector3d scene_min, scene_max;
//         for (int i = 0; i < scene_bbox.size(); i++) {
//             assert(!empty(scene_bbox(i)));
//             scene_min(i) = scene_bbox(i).lower();
//             scene_max(i) = scene_bbox(i).upper();
//         }
//         if ((scene_min.array() >= conservative_min.array()).all()
//             && (scene_max.array() <= conservative_max.array()).all()) {
//             fmt::print("========================\n");
//         }
//         std::cout << logger::fmt_eigen((scene_min - conservative_min).eval())
//                   << std::endl;
//         std::cout << logger::fmt_eigen((conservative_max - scene_max).eval())
//                   << std::endl;
//         fmt::print(
//             "{:d},{:g},{:g}\n", n, (scene_max - scene_min).norm(),
//             timer.getElapsedTime());
//     }
// }
