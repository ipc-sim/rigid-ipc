#include <logger.hpp>
#include <profiler.hpp>

#include <CLI/CLI.hpp>
#include <nlohmann/json.hpp>

#include <ccd/collision_detection.hpp>
#include <geometry/distance.hpp>
#include <io/serialize_json.hpp>

int main(int argc, char* argv[])
{
    spdlog::set_level(spdlog::level::info);

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

        ccd::EdgeVertexCandidates ev_candidates;
        ccd::detect_edge_vertex_collision_candidates(
            vertices, displacements, edges, group_ids, ev_candidates,
            ccd::DetectionMethod::BRUTE_FORCE, 0.0);
        if (ev_candidates.size() == 0) {
            csv << fmt::format("{},\n", i);
            continue;
        }
        double min_distance = -1;
        for (size_t j = 0; j < ev_candidates.size(); j++) {
            const ccd::EdgeVertexCandidate& ev_candidate = ev_candidates[j];

            double distance = ccd::geometry::point_segment_distance<double>(
                vertices.row(ev_candidate.vertex_index),
                vertices.row(edges(ev_candidate.edge_index, 0)),
                vertices.row(edges(ev_candidate.edge_index, 1)));

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
