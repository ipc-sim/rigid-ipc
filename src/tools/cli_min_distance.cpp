#include <CLI/CLI.hpp>
#include <nlohmann/json.hpp>

#include <ipc/ipc.hpp>

#include <io/read_rb_scene.hpp>
#include <io/serialize_json.hpp>
#include <physics/rigid_body_assembler.hpp>

#include <logger.hpp>

int main(int argc, char* argv[])
{
    using namespace ccd;
    using namespace physics;
    logger::set_level(spdlog::level::info);

    struct {
        std::string input_json = "";
        std::string output_csv = "";
        std::string header = "min_distance";
    } args;

    CLI::App app {
        "compute minimum distance for each step of simulation results"
    };

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
        return 1;
    }

    auto state_sequence =
        scene["animation"]["state_sequence"].get<std::vector<nlohmann::json>>();
    if (state_sequence.size() == 0) {
        return 0;
    }

    std::vector<physics::RigidBody> rbs;
    io::read_rb_scene(scene["args"]["rigid_body_problem"], rbs);
    RigidBodyAssembler bodies;
    bodies.init(rbs);

    std::stringstream csv;
    csv << fmt::format("it, {}\n", args.header);
    for (size_t i = 0; i < state_sequence.size(); ++i) {
        Poses<double> poses(bodies.num_bodies());
        assert(state_sequence[i]["rigid_bodies"].size() == bodies.num_bodies());
        for (int j = 0; j < bodies.num_bodies(); j++) {
            const auto& jrb = state_sequence[i]["rigid_bodies"][j];
            io::from_json(jrb["position"], poses[j].position);
            io::from_json(jrb["rotation"], poses[j].rotation);
        }

        Eigen::MatrixXd V = bodies.world_vertices(poses);
        ipc::Constraints constraint_set;
        ipc::construct_constraint_set(
            /*V_rest=*/V, V, bodies.m_edges, bodies.m_faces,
            /*dhat=*/1, constraint_set,
            /*ignore_internal_vertices=*/false,
            /*vertex_group_ids=*/bodies.group_ids(), bodies.m_faces_to_edges);
        double min_distance = sqrt(ipc::compute_minimum_distance(
            V, bodies.m_edges, bodies.m_faces, constraint_set));

        csv << fmt::format("{},{:.18e}\n", i, min_distance);
    }
    std::ofstream myfile;
    myfile.open(args.output_csv);
    myfile << csv.str();
    myfile.close();
}
