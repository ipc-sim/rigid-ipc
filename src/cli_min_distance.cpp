#include <logger.hpp>
#include <profiler.hpp>

#include <CLI/CLI.hpp>
#include <nlohmann/json.hpp>

#include <ipc/distance/point_edge.hpp>

#include <ccd/collision_detection.hpp>
#include <io/serialize_json.hpp>
#include <logger.hpp>

int main(int argc, char* argv[])
{
    ccd::logger::set_level(spdlog::level::info);

    struct {
        std::string input_json = "";
        std::string output_csv = "";
        std::string header = "min_distance";
    } args;

    CLI::App app { "run headless simulation" };

    app.add_option(
           "input_json,-i,--input", args.input_json,
           "JSON file with input simulation.")
        ->required();

    app.add_option(
           "output_csv,-o,--outputh", args.output_csv, "CSV file for output.")
        ->required();
    app.add_option("--header", args.header, "use as name for header.");

    try {
        app.parse(argc, argv);
    } catch (const CLI::ParseError& e) {
        return app.exit(e);
    }

    using nlohmann::json;
    std::ifstream input(args.input_json);
    json scene = json::parse(input, nullptr, false);
    if (scene.is_discarded()) {
        spdlog::error("Invalid Json file");
        exit(1);
    }

    Eigen::MatrixXi edges;
    ccd::io::from_json<int>(scene["animation"]["edges"], edges);
    Eigen::VectorXi group_ids;
    ccd::io::from_json<int>(scene["animation"]["group_id"], group_ids);

    auto& vtx_sequence = scene["animation"]["vertices_sequence"];

    std::stringstream csv;
    csv << fmt::format("it, {}\n", args.header);
    for (size_t i = 0; i < vtx_sequence.size(); ++i) {
        auto& jv = vtx_sequence[i];
        Eigen::MatrixXd vertices, displacements;
        ccd::io::from_json<double>(jv, vertices);

        assert(vertices.rows() == group_ids.rows());

        displacements.resizeLike(vertices);
        displacements.setZero();

        ipc::Candidates candidates;
        ccd::detect_collision_candidates(
            vertices, displacements, edges, Eigen::MatrixXi(0, 3), group_ids,
            ccd::CollisionType::EDGE_VERTEX, candidates,
            ccd::DetectionMethod::BRUTE_FORCE, 0.0);
        if (candidates.ev_candidates.size() == 0) {
            csv << fmt::format("{},\n", i);
            continue;
        }
        double min_distance = -1;
        for (size_t j = 0; j < candidates.ev_candidates.size(); j++) {
            const ipc::EdgeVertexCandidate& ev_candidate =
                candidates.ev_candidates[j];

            Eigen::VectorXd p = vertices.row(ev_candidate.vertex_index);
            Eigen::VectorXd s0 =
                vertices.row(edges(ev_candidate.edge_index, 0));
            Eigen::VectorXd s1 =
                vertices.row(edges(ev_candidate.edge_index, 1));
            double distance = ipc::point_edge_distance(p, s0, s1);

            if (min_distance < 0 || distance < min_distance) {
                min_distance = distance;
            }
        }
        csv << fmt::format("{},{:.18e}\n", i, min_distance);
    }
    std::ofstream myfile;
    myfile.open(args.output_csv);
    myfile << csv.str();
    myfile.close();
}
