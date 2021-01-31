#include <CLI/CLI.hpp>

#include <ccd/ccd.hpp>
#include <io/read_rb_scene.hpp>
#include <io/serialize_json.hpp>

#include <logger.hpp>

int main(int argc, char* argv[])
{
    using namespace ipc::rigid;
    set_logger_level(spdlog::level::info);

    CLI::App app { "check for collisions over simulation results" };

    std::string input_filename = "";
    app.add_option(
           "input_filename,-i,--input", input_filename,
           "JSON file with input simulation.")
        ->required();

    try {
        app.parse(argc, argv);
    } catch (const CLI::ParseError& e) {
        return app.exit(e);
    }

    using nlohmann::json;
    std::ifstream input(input_filename);
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

    int collision_types = bodies.dim() == 2
        ? CollisionType::EDGE_VERTEX
        : (CollisionType::EDGE_EDGE | CollisionType::FACE_VERTEX);

    PosesD poses_t0(bodies.num_bodies());
    assert(state_sequence[0]["rigid_bodies"].size() == bodies.num_bodies());
    for (int i = 0; i < bodies.num_bodies(); i++) {
        const auto& jrb = state_sequence[0]["rigid_bodies"][i];
        from_json(jrb["position"], poses_t0[i].position);
        from_json(jrb["rotation"], poses_t0[i].rotation);
    }

    int collision_steps = 0;
    for (size_t i = 1; i < state_sequence.size(); ++i) {
        PosesD poses_t1(bodies.num_bodies());
        assert(state_sequence[i]["rigid_bodies"].size() == bodies.num_bodies());
        for (int j = 0; j < bodies.num_bodies(); j++) {
            const auto& jrb = state_sequence[i]["rigid_bodies"][j];
            from_json(jrb["position"], poses_t1[j].position);
            from_json(jrb["rotation"], poses_t1[j].rotation);
        }

        Impacts impacts;
        detect_collisions(
            bodies, poses_t0, poses_t1, collision_types, impacts,
            DetectionMethod::HASH_GRID, TrajectoryType::RIGID);
        if (impacts.size() != 0) {
            fmt::print("{}: step {:d} failed\n", input_filename, i);
            collision_steps++;
        }

        poses_t0 = poses_t1;
    }
    fmt::print("collision_steps = {:d}\n", collision_steps);
}
