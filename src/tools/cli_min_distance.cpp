#include <CLI/CLI.hpp>
#include <nlohmann/json.hpp>

#include <ipc/ipc.hpp>

#include <io/read_rb_scene.hpp>
#include <io/serialize_json.hpp>
#include <physics/rigid_body_assembler.hpp>

#include <logger.hpp>

int main(int argc, char* argv[])
{
    using namespace ipc;
    using namespace ipc::rigid;
    set_logger_level(spdlog::level::info);

    CLI::App app {
        "compute minimum distance for each step of simulation results"
    };

    std::string input_json = "";
    app.add_option(
           "input_json,-i,--input", input_json,
           "JSON file with input simulation.")
        ->required();

    std::string output_csv = "";
    app.add_option(
           "output_csv,-o,--outputh", output_csv, "CSV file for output.")
        ->required();

    std::string header = "min_distance";
    app.add_option("--header", header, "use as name for header.");

    try {
        app.parse(argc, argv);
    } catch (const CLI::ParseError& e) {
        return app.exit(e);
    }

    using nlohmann::json;
    std::ifstream input(input_json);
    json scene = json::parse(input, nullptr, false);
    if (scene.is_discarded()) {
        spdlog::error("Invalid Json file");
        return 1;
    }

    auto state_sequence =
        scene["animation"]["state_sequence"].get<std::vector<nlohmann::json>>();
    if (state_sequence.size() == 0) {
        return 0;
    }

    std::vector<RigidBody> rbs;
    read_rb_scene(scene["args"]["rigid_body_problem"], rbs);
    RigidBodyAssembler bodies;
    bodies.init(rbs);

    std::stringstream csv;
    csv << fmt::format("it, {}\n", header);
    for (size_t i = 0; i < state_sequence.size(); ++i) {
        PosesD poses(bodies.num_bodies());
        assert(state_sequence[i]["rigid_bodies"].size() == bodies.num_bodies());
        for (int j = 0; j < bodies.num_bodies(); j++) {
            const auto& jrb = state_sequence[i]["rigid_bodies"][j];
            from_json(jrb["position"], poses[j].position);
            from_json(jrb["rotation"], poses[j].rotation);
        }

        Eigen::MatrixXd V = bodies.world_vertices(poses);

        const Eigen::VectorXi& group_ids = bodies.group_ids();
        auto can_collide = [&group_ids](size_t vi, size_t vj) {
            return group_ids[vi] != group_ids[vj];
        };

        Constraints constraint_set;
        construct_constraint_set(
            /*V_rest=*/V, V, bodies.m_edges, bodies.m_faces,
            /*dhat=*/1, constraint_set, bodies.m_faces_to_edges, /*dmin=*/0,
            BroadPhaseMethod::HASH_GRID,
            /*ignore_internal_vertices=*/false, can_collide);
        double min_distance = sqrt(compute_minimum_distance(
            V, bodies.m_edges, bodies.m_faces, constraint_set));

        csv << fmt::format("{},{:.18e}\n", i, min_distance);
    }
    std::ofstream myfile;
    myfile.open(output_csv);
    myfile << csv.str();
    myfile.close();
}
