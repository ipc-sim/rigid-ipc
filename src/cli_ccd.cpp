#include <logger.hpp>
#include <profiler.hpp>

#include <CLI/CLI.hpp>

#include <ccd/collision_detection.hpp>
#include <io/serialize_json.hpp>
#include <opt/distance_barrier_constraint.hpp>

int main(int argc, char* argv[])
{
    spdlog::set_level(spdlog::level::info);

    struct {
        std::string input_json = "";
        std::string output_csv = "";
    } args;

    CLI::App app { "run headless simulation" };

    app.add_option("input_json,-i,--input", args.input_json,
           "JSON file with input simulation.")
        ->required();

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
    int collision_steps = 0;
    for (size_t i = 0; i < vtx_sequence.size() - 1; ++i) {
        auto& jv = vtx_sequence[i];
        auto& jv2 = vtx_sequence[i + 1];
        Eigen::MatrixXd vertices, displacements;
        ccd::io::from_json<double>(jv, vertices);
        ccd::io::from_json<double>(jv2, displacements);
        displacements = displacements - vertices;

        assert(vertices.rows() == group_ids.rows());

        ccd::EdgeVertexImpacts ev_impacts;
        ccd::detect_edge_vertex_collisions(vertices, displacements, edges,
            group_ids, ev_impacts, ccd::DetectionMethod::BRUTE_FORCE);
        if (ev_impacts.size() != 0) {
            std::cout << args.input_json << ": step " << i << " failed"
                      << std::endl;
            collision_steps+=1;
        }

    }
    std::cout << "collision_steps = " << collision_steps << std::endl;
}
